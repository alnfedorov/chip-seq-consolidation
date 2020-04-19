import os
import logging
from logging import Handler, LogRecord, Formatter


def config_logging(logfolder: str, level: str = "NOTSET"):
    assert os.path.isdir(logfolder)
    logging.basicConfig(level=level, format="%(levelname)s: %(name)s: %(message)s")
    logger = logging.getLogger()

    handler = PerFileHandler(logfolder)
    handler.setLevel(level)
    handler.setFormatter(Formatter("%(asctime)s: %(levelname)s: %(message)s"))
    logger.addHandler(handler)


class PerFileHandler(Handler):
    def __init__(self, logfolder, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.logfolder = logfolder
        self.__files = {}

    def emit(self, record: LogRecord) -> None:
        file = os.path.join(self.logfolder, record.name)
        if file not in self.__files:
            self.__files[file] = open(file, 'a')
        print(self.format(record), file=self.__files[file], flush=True)

    def flush(self) -> None:
        for f in self.__files.values():
            f.close()
        self.__files.clear()
