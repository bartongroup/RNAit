#!/usr/bin/env python

import cgi
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio.Application import ApplicationError
from tempfile import mkdtemp
import io
import re
import primer3
import pprint
import shutil
import textwrap
import os
from jinja2 import Environment, FileSystemLoader, select_autoescape
import yaml

import cgitb
cgitb.enable(format='text')


def application(environ, start_response):
    start_response('200 OK', [('Content-Type', 'text/html')])

    config_file = (
        os.path.dirname(
            os.path.realpath(__file__)) +
        '/../etc/RNAit.yaml')
    with open(config_file) as s:
        config = yaml.safe_load(s)

    RNAit_dir = config.get('RNAit_dir')
    db_dir = config.get('db_dir')
    tmp_dir = config.get('tmp_dir')

    params = get_params(environ)

    if ('error' in params):
        return(get_error_page(RNAit_dir, params.get('error'), 'submission'))

    seq = params.get('seq')
    db = ('%s%s' % (db_dir, params.get('database')))
    primers, error = get_primer_pairs(params)
    # capture query parameters for display on results page
    query_info = {
        'query_seq': seq.id,
        'query_length': len(seq.seq),
        'melting_temp': params.get('melting_temp'),
        'product_size': params.get('product_size'),
        'database': params.get('database'),
        'stringency': "%s - %s" % (params.get('string_min'), params.get('string_max')),
        'subunit_length': params.get('subunit_length'),
    }

    if error:
        return(get_error_page(RNAit_dir, error, 'runtime'))

    if (len(primers) == 0):
        return(get_error_page(RNAit_dir, 'No suitable primers found', 'runtime'))
        encode = html.encode('UTF-8')
        return(encode)

    blast_results = []
    string_min = int(params.get('string_min'))
    string_max = int(params.get('string_max'))
    subunit_length = int(params.get('subunit_length'))

    for pair in primers:
        product = get_pcr_product(seq, pair)
        blast_result, error = blast_product(
            product, tmp_dir, db, string_min, string_max, subunit_length)
        if (error):
            return(get_error_page(RNAit_dir, error, 'runtime'))
        blast_results.append(blast_result)
    html = get_output_page(query_info, primers, RNAit_dir, blast_results)

    return [html]

# get_params
#
# parses form parameters following form submission
# Sequences entered via the 'seqpaste' field are parsed
# as SeqRecord objects via SeqIO
#
# required args: environ - environment dictionary
#
# returns: params - dictionary of parsed parameters


def get_params(environ):

    post_env = environ.copy()
    post_env['QUERY_STRING'] = ''
    post = cgi.FieldStorage(
        fp=environ['wsgi.input'],
        environ=post_env,
        keep_blank_values=True)

    params = {}
    for f in post.list:
        if f.name == 'seqpaste' and f.value != '':
            seqH = io.StringIO(f.value)
            try:
                record = SeqIO.read(seqH, 'fasta')
                params['seq'] = record
            except ValueError:
                params['error'] = 'the entered sequence does not appear to be valid fasta format'
        elif f.name == 'upload' and post.getvalue('seqpaste') == '':
            raw_fasta = f.value.decode("utf-8")
            seqH = io.StringIO(raw_fasta)
            try:
                record = SeqIO.read(seqH, 'fasta')
                params['seq'] = record
            except ValueError:
                params['error'] = 'the uploaded sequence does not appear to be valid fasta format'
        else:
            params[f.name] = f.value

    return(params)

# get_primer_pairs
#
# runs primer3 using primer3-py bindings against the provided sequence according to the requested
# melting temperature and product size parameters
#
# required args: params - dictionary of parsed form parameters
#
# returns: primers - list of primer pair dictionaries
#          error - runtime error (string)


def get_primer_pairs(params):

    seq_args = {
        'SEQUENCE_ID': params.get('seq').id,
        'SEQUENCE_TEMPLATE': str(params.get('seq').seq),
    }

    global_args = {
        'PRIMER_TASK': 'generic',
        'PRIMER_PRODUCT_SIZE_RANGE': [[int(params.get('product_min')), int(params.get('product_max'))]],
        'PRIMER_OPT_TM': int(params.get('melting_temp')),
    }

    primers = {}
    try:
        primers = primer3.bindings.designPrimers(seq_args, global_args)
    except OSError as error:
        return(primers, error)

    pair_count = int(primers.get('PRIMER_PAIR_NUM_RETURNED'))
    pairs = []
    for i in range(pair_count):
        formatted_product = get_formatted_product(
            str(params.get('seq').seq), primers, i)

        pair = {
            'LEFT_START': (primers.get('PRIMER_LEFT_' + str(i)))[0],
            'LEFT_LENGTH': (primers.get('PRIMER_LEFT_' + str(i)))[1],
            'RIGHT_START': (primers.get('PRIMER_RIGHT_' + str(i)))[0],
            'RIGHT_LENGTH': (primers.get('PRIMER_RIGHT_' + str(i)))[1],
            'LEFT_SEQ': primers.get('PRIMER_LEFT_' + str(i) + '_SEQUENCE'),
            'RIGHT_SEQ': primers.get('PRIMER_RIGHT_' + str(i) + '_SEQUENCE'),
            'LEFT_GC': primers.get('PRIMER_LEFT_' + str(i) + '_GC_PERCENT'),
            'RIGHT_GC': primers.get('PRIMER_RIGHT_' + str(i) + '_GC_PERCENT'),
            'LEFT_MELTING': "%.2f" % primers.get('PRIMER_LEFT_' + str(i) + '_TM'),
            'RIGHT_MELTING': "%.2f" % primers.get('PRIMER_RIGHT_' + str(i) + '_TM'),
            'LEFT_END_STAB': primers.get('PRIMER_LEFT_' + str(i) + '_END_STABILITY'),
            'RIGHT_END_STAB': primers.get('PRIMER_RIGHT_' + str(i) + '_END_STABILITY'),
            'PRODUCT_SIZE': primers.get('PRIMER_PAIR_' + str(i) + '_PRODUCT_SIZE'),
            'COMP_END': primers.get('PRIMER_PAIR_' + str(i) + '_COMPL_END'),
            'PRODUCT': formatted_product,
        }
        pairs.append(pair)

    return(pairs, None)

# get_formatted_product
#
# Adds HTML highlighting to region of sequence to be amplified by primers
#
# required args: seq - sequence
#                primers - dictionary of primers returned from primer3.binding.designPrimers
#                i - index of primer pair to use
#
# returns: formatted_seq - sequence with html formatting applied


def get_formatted_product(seq, primers, i):

    start = (primers.get('PRIMER_LEFT_' + str(i)))[0]
    end = (primers.get('PRIMER_RIGHT_' + str(i)))[0]

    lines = textwrap.wrap(seq, width=60)
    count = 0
    formatted_seq = ''
    for line in lines:
        max_length = count + 60
        line_start = count
        count = count + len(line)
        num_spaces = max_length - count + 4
        spaces = '&nbsp;' * num_spaces

        # add a <span> at the beginning of the product
        if(start > line_start and start < count):
            offset = start - line_start
            left_flank = line[:offset]
            product = line[offset:]
            line = "%s%s%s%s" % (
                left_flank, '<span style="color:red">', product, '</span>')

        # wrap lines completely within product in <span>s
        if(end > line_start and start < line_start and end > count):
            line = "%s%s%s" % ('<span style="color:red">', line, '</span>')

        # add a <span> around the end of the product
        if(end > line_start and end <= count):
            offset = end - line_start
            product = line[:offset + 1]
            right_flank = line[offset + 1:]
            line = "%s%s%s%s" % ('<span style="color:red">',
                                 product, '</span>', right_flank)

        formatted_seq = "%s%s%s%s<br/>" % (formatted_seq, line, spaces, count)

    return(formatted_seq)

# get_pcr_product
#
# Isolates the subsequence represnting the pcr product for a primer pair
#
# required args: seq - Bio.seqRecord object
#                pair - dictionary of primer pair info
#
# returns: product -  Bio.seqRecord object


def get_pcr_product(seq, pair):

    # both primer3 and biopython use 0 based co-ordinates so we should be good
    # to go....
    start = pair.get('LEFT_START')
    end = pair.get('RIGHT_START') + 1
    product = seq[start:end]
    return(product)

# blast_product
#
# Blasts pcr product against organism genome database to identify
#
# required args: product - Bio:seqRecord object representing pcr product
#                db - blast database name
#
# returns: blast_data - dictionary containing 'record' (blast_record
# object), alignment status etc.


def blast_product(product, tmp_dir, db, string_min,
                  string_max, subunit_length):
    blast_dir = mkdtemp(dir=tmp_dir)
    queryFileName = blast_dir + '/query'
    outFileName = blast_dir + '/output.xml'

    SeqIO.write(product, queryFileName, 'fasta')
    cline = NcbiblastnCommandline(
        cmd='blastn',
        query=queryFileName,
        out=outFileName,
        outfmt=5,
        db=db,
        evalue=0.01)
    stderr = ''
    try:
        stdout, stderr = cline()
    except ApplicationError as err:
        return('', err.stderr)

    result_handle = open(outFileName)
    status = ''

    blast_record = NCBIXML.read(result_handle)
    midline_regex = re.compile(r"\|{20,}")
    alignment_status = ''
    self_alignments = []
    conflicting_alignments = []
    matching_alignments = []
    reasons = []

    # flag for tracking conflicting hits identified
    conflicting = 0
    # flag for tracking hits with match exceeding subunit length
    matching = 0

    for alignment in blast_record.alignments:

        alignment_data = {
            'accession': alignment.hit_id,
            'description': alignment.hit_def,
            'subj_length': alignment.length,
        }
        hsp_count = 0
        hsp_idents = []
        # lengths of consecutive bases...
        hsp_match_lengths = []
        hsp_alignments = []
        hsp_hit_lengths = []

        # Original RNAit implementation reports single value for identity, which
        # is tricky without tiling HSPs We'll use some slightly different critera
        # here

        # >1 hsp suggests the alignment is to a repetitive sequence which is
        # unlikely to amplify cleanly so mark these as bad

        for hsp in alignment.hsps:
            hsp_count += 1
            # check for matches of >20bp identity by checking for stretches of
            # >20 '|' characters in the HSP midline
            match = midline_regex.search(hsp.match)
            match_len = match.end() - match.start()
            ident = (hsp.identities / hsp.align_length)
            hsp_idents.append(ident)
            hsp_match_lengths.append(match_len)
            hsp_hit_lengths.append(hsp.align_length)

            # pretty format alignment
            text_alignment = format_alignment(hsp)
            hsp_alignments.append(text_alignment)

        if (hsp_count == 1):
            length_cov = hsp_match_lengths[0] / blast_record.query_letters
            if (hsp_idents[0] > 0.99 and length_cov == 1):
                alignment_status = 'Self alignment'
                self_alignments.append(alignment_data)
            elif (hsp_idents[0] > string_min and hsp_idents[0] < string_max):
                alignment_status = 'Conflicting hits'
                conflicting_alignments.append(alignment_data)
                conflicting += 1
                reasons.append("Identity is %s" % (hsp_idents[0]))
            elif (hsp_match_lengths[0] > subunit_length):
                alignment_status = 'Match exceeding subunit length'
                matching_alignments.append(alignment_data)
                matching += 1
                reasons.append(
                    "%s bp identical sequence" %
                    (hsp_match_lengths[0]))
            else:
                alignment_status = 'Good'
        else:
            alignment_status = 'Multiple HPSs'

        hsp_idents = list(map(format_ident, hsp_idents))
        alignment_data['status'] = alignment_status
        alignment_data['reasons'] = reasons
        alignment_data['hsps'] = hsp_count
        alignment_data['ident'] = ";".join(map(str, hsp_idents))
        alignment_data['hsp_alignments'] = hsp_alignments
        alignment_data['hsp_hit_lengths'] = ";".join(map(str, hsp_hit_lengths))

    if conflicting:
        primer_status = 'Bad'
    elif matching:
        primer_status = 'Bad'
    else:
        primer_status = 'Suitable'

    blast_data = {
        'record': blast_record,
        'primer_status': primer_status,
        'self_alignments': self_alignments,
        'conflicting_alignments': conflicting_alignments,
        'matching_alignments': matching_alignments,
    }
    shutil.rmtree(blast_dir)

    return(blast_data, None)

# get_output_page
#
# Generates HTML output based on jinja template
#
# required args: query_info - dictionary of query data
#                primers - dictionary of primers produced by get_primer_pairs
#                blast_data - dictionary of blast results generated by blast_product
#
# returns: page - HTML page


def get_output_page(query_info, primers, RNAit_dir, blast_results):
    env = Environment(
        loader=FileSystemLoader(RNAit_dir + '/templates'), autoescape=select_autoescape(['html', 'xml'])
    )
    template = env.get_template('result_page.html')
    html = template.render(
        query_info=query_info,
        primers=primers,
        blast=blast_results)
    encode = html.encode('UTF-8')
    return(encode)

# get_error_page
#
# Generates HTML error page based on jinga template
#
# required args: RNAit_dir - path to RNAit installation (string)
#                error - Error message to render (string)
#                type - submission or runtime (string)
#                       'submission' error warns regaring invalid inputs


def get_error_page(RNAit_dir, error, type):
    env = Environment(
        loader=FileSystemLoader(RNAit_dir + '/templates'), autoescape=select_autoescape(['html', 'xml'])
    )
    template = env.get_template('error_page.html')
    html = template.render(error=error, type=type)
    encode = html.encode('UTF-8')
    return(encode)

# format_ident
#
# Converts proportion of identities (i.e. 0.93) to a percentage
#
# required args: val (str)
#
# returns: percent (float)


def format_ident(val):
    percent = ("%.2f" % (float(val) * 100))
    return(percent)

# format_alignment
#
# Produces a text alignment from hsp with HTML linkbreaks/spaces
# for rendering in a <pre>
#
# requred args: hsp - Bio.Blast.Record.HSP
#
# returns: alignment - string


def format_alignment(hsp):
    query_start = hsp.query_start
    sbjct_start = hsp.sbjct_start
    offset = 0
    alignment_lines = []

    line = "Score: %.2f; bits: %.2f; e-value: %.2f\n" % (
        hsp.score, hsp.bits, hsp.expect)
    alignment_lines.append(line)
    alignment_lines.append("\n")

    for i in range(0, hsp.align_length, 75):

        qline = ("Query: %s %s" %
                 (str(query_start + offset).rjust(4), hsp.query[i:i + 75]))
        midline = ("            %s" % (hsp.match[i:i + 75]))
        hline = ("Sbjct: %s %s" %
                 (str(sbjct_start + offset).rjust(4), hsp.sbjct[i:i + 75]))

        alignment_lines.append(qline)
        alignment_lines.append(midline)
        alignment_lines.append(hline)
        alignment_lines.append('')
        offset += 75

    alignment = "<br/>".join(alignment_lines)
    alignment = alignment.replace(' ', '&nbsp;')
    return(alignment)
