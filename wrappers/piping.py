import tempfile
import asyncio
import os
import fcntl
from typing import Awaitable, Callable, Collection
import collections
from inspect import Signature, BoundArguments
from functools import partialmethod


class PipingAwaitable(Awaitable[str]):
    def __init__(self, func: Callable, bargs: BoundArguments, writearg: str, enforce_file: bool,
                 pawaitables: {str: 'PipingAwaitable'}):
        self.underlying = func
        self.bargs = bargs
        self.awaited = False

        assert writearg in bargs.arguments.keys()
        assert all(k in bargs.arguments.keys() for k in pawaitables)

        self.writearg = writearg
        self.pawaitables = pawaitables

        self.mode = "file"
        self.enforce_file = enforce_file

        # Always respect user passed paths
        if self.bargs.arguments[self.writearg] is not None:
            self.enforce_file = True
            self.consumers = float("inf")
        else:
            self.consumers = 0

        self.result = None
        self.status = "initialized"
        self.coro = None

    def ensure_file(self):
        self.consumers += 1
        if self.status == "finished" and self.mode == "file" and os.path.isfile(self.result):
            return

        assert self.result is None
        self.mode = "file"
        self.enforce_file = True

    def maybe_pipe(self):
        self.consumers += 1

        if self.status == "finished" and self.mode == "file" and os.path.isfile(self.result):
            return

        assert self.result is None
        if self.consumers == 1 and not self.enforce_file:
            self.mode = "pipe"
        else:
            self.mode = "file"
            self.enforce_file = True

    async def launch(self):
        # Intermediates with file output can be used multiple times
        if self.mode == "file" and self.status != "initialized":
            return self.result

        assert self.status == "initialized"
        self.status = "launched"
        # 1. Substitute paths with real files
        for k, pipe in self.pawaitables.items():
            args = await pipe.launch()
            self.bargs.arguments[k] = args

        # 2. Create file and launch process
        if self.bargs.arguments[self.writearg] is None:
            self.bargs.arguments[self.writearg] = tempfile.mkstemp()[1]
        self.result = self.bargs.arguments[self.writearg]
        if self.mode == "file":
            # file must exist before downstream program start
            assert await self.underlying(*self.bargs.args, **self.bargs.kwargs) == self.result
        else:
            assert self.mode == "pipe"
            os.remove(self.result)  # file shouldn't exist before fifo call
            os.mkfifo(self.result)
            file = os.open(self.result, os.O_RDONLY | os.O_NONBLOCK)
            fcntl.fcntl(file, 1031, 200 * 1024 * 1024 // 4096)  # F_SETPIPE_SZ(1031) to 200MB
            os.close(file)
            self.coro = asyncio.create_task(self.underlying(*self.bargs.args, **self.bargs.kwargs))
        return self.result

    async def finalize(self):
        # Intermediates with file output can be used multiple times
        if self.status == "finished" and self.mode == "file":
            return

        assert self.result is not None
        self.status = "finalized"
        # 1. Await finishing of all ancestors
        for pipe in self.pawaitables.values():
            await pipe.finalize()

        # 2. Await launched process
        if self.mode == "pipe":
            assert await self.coro == self.result

    async def cleanup(self, intermediate: bool = True):
        assert self.result is not None
        if self.status != "finished":
            assert self.status == "finalized"
            # Await cleanup of all ancestors
            for pipes in self.pawaitables.values():
                await pipes.cleanup(intermediate=True)
            self.consumers -= 1
            self.status = "finished"

        # 2. Remove pipe or tmp file if possible
        if intermediate and self.consumers <= 0:
            os.remove(self.result)
        return self.result

    def __await__(self):
        """If this method is called, then it is the end of the chain. And user needs its output"""
        assert self.result is None, \
            "PipingAwaitable is being used twice. It is prohibited to split pipeline and reuse intermediates " \
            "when pipeline finishes. Explicitly await intermediate results if needed."

        # If this file is called with await -> we must save the result:
        self.mode = 'file'
        self.consumers = float("inf")

        # 1. Launch processing
        yield from self.launch().__await__()
        # 2. Await finishing
        yield from self.finalize().__await__()
        # 3. Do cleanup
        return (yield from self.cleanup(intermediate=False).__await__())


class CompositePipingAwaitable:
    def __init__(self, wrapped: Collection[PipingAwaitable]):
        self.mode = CompositePipingAwaitable.mode(wrapped)
        assert self.mode
        self.wrapped = wrapped

    def _delegate(self, method, **kwargs):
        if self.mode == "sequence":
            return [getattr(v, method)(**kwargs) for v in self.wrapped]
        elif self.mode == "set":
            return {getattr(v, method)(**kwargs) for v in self.wrapped}
        elif self.mode == "map":
            return {k: getattr(v, method)(**kwargs) for k, v in self.wrapped.items()}

    async def _adelegate(self, method, **kwargs):
        if self.mode == "sequence":
            values = await asyncio.gather(*[getattr(v, method)(**kwargs) for v in self.wrapped])
            return list(values)
        elif self.mode == "set":
            values = await asyncio.gather(*[getattr(v, method)(**kwargs) for v in self.wrapped])
            return set(values)
        elif self.mode == "map":
            keys, values = zip(*self.wrapped.items())
            values = await asyncio.gather(*[getattr(v, method)(**kwargs) for v in self.wrapped])
            return {k: v for k, v in zip(keys, values)}

    ensure_file = partialmethod(_delegate, method="ensure_file")
    maybe_pipe = partialmethod(_delegate, method="maybe_pipe")
    launch = partialmethod(_adelegate, method="launch")
    finalize = partialmethod(_adelegate, method="finalize")
    cleanup = partialmethod(_adelegate, method="cleanup")

    @staticmethod
    def mode(awaitables):
        if isinstance(awaitables, collections.abc.Sequence) and all(isinstance(x, PipingAwaitable) for x in awaitables):
            return "sequence"
        elif isinstance(awaitables, collections.abc.Set) and all(isinstance(x, PipingAwaitable) for x in awaitables):
            return "set"
        elif isinstance(awaitables, collections.abc.Mapping) and \
                all(isinstance(x, PipingAwaitable) for x in awaitables.values()):
            return "map"
        return None


def pipe(writearg: str, **pipe_kwargs):
    """
    Если какие-то из аргументов в pipe_kwargs можно сделать трубой -
    сигнализирует об этом соответствующим PipingAwaitable.
    """
    writearg, writetomode = writearg
    assert all(k in ("p", "f") for k in pipe_kwargs.values()) and writetomode in ("p", "f")

    def wrapper(func):
        sig: Signature = Signature.from_callable(func, follow_wrapped=True)

        # Check prerequisites
        assert sig.return_annotation == str, \
            "piping works only with functions that returns a single arguments"
        assert all(k in sig.parameters.keys() for k in pipe_kwargs), \
            f"Failed to match pipe-enabled keys with function signature " \
            f"{list(sig.parameters.keys())} != {list(pipe_kwargs.keys())}"
        assert writearg in sig.parameters.keys(), \
            f"Missing writeto {writearg} in the function parameters list {list(sig.parameters.keys())}"

        def newfunc(*args, **kwargs):
            barguments = sig.bind_partial(*args, **kwargs)
            barguments.apply_defaults()
            allargs = barguments.arguments
            # Check if we can try to make pipe from pipe-enabled args
            piping_awaitables = {}
            for key, mode in pipe_kwargs.items():
                if CompositePipingAwaitable.mode(allargs[key]) is not None:
                    value = CompositePipingAwaitable(allargs[key])
                elif isinstance(allargs[key], PipingAwaitable):
                    value = allargs[key]
                else:
                    continue
                if mode == 'f':
                    value.ensure_file()
                    piping_awaitables[key] = value
                else:
                    assert mode == 'p'
                    value.maybe_pipe()
                    piping_awaitables[key] = value
            return PipingAwaitable(func, barguments, writearg, writetomode == 'f', pawaitables=piping_awaitables)

        return newfunc

    return wrapper
