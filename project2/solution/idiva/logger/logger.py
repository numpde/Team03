# RA, 2020-12-01

import logging
import pathlib
import datetime

# Based in part on
# https://dev.to/joaomcteixeira/setting-up-python-logging-for-a-library-app-6ml
# https://github.com/numpde/transport/blob/9b53b7c/pt2pt/20181019-study1/helpers/commons.py

LOG_FILE = (pathlib.Path(__file__).parent / "logs")
LOG_FILE.mkdir(exist_ok=True, parents=True)

for f in sorted(LOG_FILE.glob("*-*-*.log"))[:-10]:
    os.remove(f)

LOG_FILE = LOG_FILE / datetime.datetime.now(tz=datetime.timezone.utc).strftime("%Z-%Y%m%d-%H%M%S")
LOG_FILE = LOG_FILE.with_suffix(".log")


class _:
    import logging.config
    logging.config.dictConfig(dict(
        version=1,
        formatters={
            'sparse': {
                'format': "[%(asctime)s] %(message)s",
                'datefmt': "%H:%M:%S %Z",
            },
            'verbose': {
                'format': "[%(levelname).1s][%(asctime)s][%(pathname)s][%(funcName)s:%(lineno)d]\n%(message)s",
                'datefmt': "%Y%m%d %H:%M:%S %Z",
            },
        },
        handlers={
            'h': {
                'class': "logging.StreamHandler",
                'formatter': "sparse",
                'level': logging.INFO,
                # Note: progressbar uses stderr
                'stream': "ext://sys.stderr",
            },
            'f': {
                'class': "logging.FileHandler",
                'formatter': "verbose",
                'level': logging.DEBUG,
                'filename': LOG_FILE,
            }
        },
        root={
            'handlers': ['h', 'f'],
            'level': logging.DEBUG,
        }
    ))


log = logging.getLogger(__name__)
