#!/usr/bin/env python
import sys
from rq import Connection, Worker
import os

from retrobiocat_web.app.app import create_app

scheduler = os.environ.get('SCHEDULER') or False
production_mode = os.environ.get('PRODUCTION') or False

if __name__ == '__main__':
    app = create_app(use_talisman=production_mode)
    app.app_context().push()

    with Connection(app.redis):
        qs = sys.argv[1:] or ['auto_jobs', 'tasks', 'network', 'pathway', 'db', 'process_blasts', 'alignment', 'blast', 'preprocess']

        w = Worker(qs, log_job_description=False)
        w.work(with_scheduler=scheduler)

