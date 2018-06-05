#!/usr/bin/env python

from cgi import parse_qs
import cgitb
cgitb.enable()

def application (environ,start_response):
    start_response('200 OK',[('Content-Type','text/html')])
    
    params = get_params(environ)
    print(params)
    
    html = b'seqpaste:'
        
    return [html]


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
        params[key.decode("utf-8")]=val
        
    #seqpaste_val = d.get(b'seqpaste',[''])[0].decode("utf-8")
    #print(seqpaste_val)f
#    melting_temp = d.get('melting_temp',[''])[0]
#    product_min = d.get('product_min',[''])[0]
 #   product_max = d.get('product_max',[''])[0]
  #  string_min = d.get('string_min',[''])[0]
   # string_max = d.get('string_max',[''])[0]
    #subunit_length = d.get('subunit_length',[''])[0]
  #  database = d.get('database',[''])[0]
    
    
    return(params)
