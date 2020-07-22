#!/usr/bin/env python
import sys
from rq import Connection, Worker

from retrobiocat_web.app.app import create_app

#retrorules = RetroRules(None)

if __name__ == '__main__':
    app = create_app()

    #retrorules.load()
    #app.retrorules_rxns = retrorules.retrorules_rxns
    #app.retrorules_db = retrorules.retrorule_db

    app.app_context().push()

    with Connection(app.redis):
        qs = sys.argv[1:] or ['tasks', 'network', 'pathway', 'db']

        w = Worker(qs, log_job_description=False)
        w.work()

