#!/usr/bin/env python

from cgi import parse_qs
from Bio import SeqIO
import io

import cgitb
cgitb.enable()

def application (environ,start_response):
    start_response('200 OK',[('Content-Type','text/html')])
    
    params = get_params(environ)
    if ('error' in params):
        html = params['error']
        
        encode=html.encode('UTF-8')
        return(encode)
    
    html=''
    encode=html.encode('UTF-8')
        
    return [encode]


# get_params
#
# parses form parameters following form submission
# Sequences entered via the 'seqpaste' field are parsed
# as SeqRecord objects via SeqIO
#
# required args: environ - environment dictionary
#
# returns: params - dictionary of parsed parameters
            
# TODO: validation of submitted values

def get_params(environ):
    
    try:
        request_body_size = int(environ.get('CONTENT_LENGTH',0))
    except (ValueError):
        request_body_size = 0
    request_body = environ['wsgi.input'].read(request_body_size) 
    d = parse_qs(request_body)
    
    params={}
    for key in d.keys():
        val = d.get(key,[''])[0].decode("utf-8")
        if (key.decode("utf-8")=='seqpaste'):
            seqH = io.StringIO(val)
            try:
                record = SeqIO.read(seqH,'fasta')
                params['seq']=record
            except ValueError:
                params['error']='The entered sequence does not appear to be valid fasta format'
        else:
            params[key.decode("utf-8")]=val
        
    return(params)
