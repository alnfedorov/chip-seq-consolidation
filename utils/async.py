import asyncio
from typing import Awaitable


async def batched_gather(tasks: [Awaitable], batch_size: int, verbose=False):
    results = []
    while tasks:
        if verbose:
            print("tasks left ", len(tasks))
        batch = tasks[:batch_size]
        if len(batch) == 0:
            break
        batch = await asyncio.gather(*batch)
        tasks = tasks[batch_size:]
        results.extend(batch)
    return results
