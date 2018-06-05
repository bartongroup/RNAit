#!/usr/bin/env bash

source activate web

/Users/jabbott/miniconda3/envs/web/bin/nginx -s stop &
/Users/jabbott/miniconda3/envs/web/bin/uwsgi --ini /Users/jabbott/Development/DAG_Website/etc/uwsgi.conf --stop /Users/jabbott/miniconda3/envs/web/var/run/uwsgi.pid &

