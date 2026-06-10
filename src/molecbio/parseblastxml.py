#!/usr/bin/python
"""
Purpose: To parse the XML output of some BLAST run, collecting related results and grouping them together in a readable format. Also allows filtering of the results based on several common metrics, and saving the filtered results to a file.
Author: Dave Curran
Version: 0.1
Was written for the output of BLAST v2.10.0.
To get the correct output from the command line version of BLAST use -outfmt 5; from the online version of BLAST select "Download All", and then "XML".
Note: hits that are common to multiple queries will be merged into a single object based on their given IDs. If the subject sequences are unnamed, this will not happen. If multiple query or subject sequences share the same ID but are actually different, they will be inappropriately combined. If multiple copies of query or subject sequences have been included, they will be ignored.
"""

# TODO:
# - Create custom errors, replace all exit() calls with them.

import os
import xml.etree.ElementTree as ET
from collections import OrderedDict

def results_from_file(file_path):
    path = os.path.realpath(file_path)
    if not os.path.isfile(path):
        print('Error: no file found at "{}"'.format(path))
        exit()
    with open(path) as f:
        res = BlastResults(f.read())
    return res
def results_from_string(xml_data):
    return BlastResults(xml_data)

class BlastResults(object):
    def __init__(self, xml_data):
        # #  Public attributes  # #
        self.blast_program = ''
        self.run_parameters = {}
        self.queries = OrderedDict()
        self.hits = OrderedDict()
        # #  Private attributes  # #
        self.root = None
        self.filter_criteria = set(["max_e_value", "min_identity", "min_positive", "min_query_coverage", "min_hit_coverage", "max_gap_coverage", "min_hsp_length", "min_bitscore", "min_score"])
        self.hsp_keys = set(["e_value", "percent_identity", "percent_positive", "query_coverage", "hit_coverage", "gap_coverage", "hsp_length", "bit_score", "score"])
        self.best_is_low = set(["e_value", "gap_coverage"])
        # #  Parse the data  # #
        self.parse_xml_data(xml_data)

    # # #  Public methods  # # #
    def get_hits_per_query(self, max_e_value=None, min_identity=None, min_positive=None, min_query_coverage=None, min_hit_coverage=None, max_gap_coverage=None, min_hsp_length=None, min_bitscore=None, min_score=None, num_matches=None, sort=None):
        # Returns a list of BlastSequence objects [(query1, OrderedDict{hit1:[hsp1, hsp2, ...], ...}), ...]. The number of hits in the OrderedDict is specified by 'num_matches', and the order by 'sort'.
        seqs = []
        for query in self.queries.values():
            matches = OrderedDict()
            hsps = query.filter_hsps(max_e_value=max_e_value, min_identity=min_identity, min_positive=min_positive, min_query_coverage=min_query_coverage, min_hit_coverage=min_hit_coverage, max_gap_coverage=max_gap_coverage, min_hsp_length=min_hsp_length, min_bitscore=min_bitscore, min_score=min_score)
            if sort != None:
                hsps.sort(key=lambda hsp: getattr(hsp, sort))
                if sort not in self.best_is_low:
                    hsps.reverse()
            for hsp in hsps:
                if num_matches != None and len(matches) == num_matches and hsp.hit not in matches:
                    break
                matches.setdefault(hsp.hit, []).append(hsp)
            if len(matches) > 0:
                seqs.append((query, matches))
        return seqs
    def get_queries_per_hit(self, max_e_value=None, min_identity=None, min_positive=None, min_query_coverage=None, min_hit_coverage=None, max_gap_coverage=None, min_hsp_length=None, min_bitscore=None, min_score=None, num_matches=None, sort=None):
        # Returns a list of BlastSequence objects [(hit1, OrderedDict{query1:[hsp1, hsp2, ...], ...}), ...]. The number of queries in the OrderedDict is specified by 'num_matches', and the order by 'sort'.
        seqs = []
        for hit in self.hits.values():
            matches = OrderedDict()
            hsps = hit.filter_hsps(max_e_value=max_e_value, min_identity=min_identity, min_positive=min_positive, min_query_coverage=min_query_coverage, min_hit_coverage=min_hit_coverage, max_gap_coverage=max_gap_coverage, min_hsp_length=min_hsp_length, min_bitscore=min_bitscore, min_score=min_score)
            if sort != None:
                hsps.sort(key=lambda hsp: getattr(hsp, sort))
                if sort not in self.best_is_low:
                    hsps.reverse()
            for hsp in hsps:
                if num_matches != None and len(matches) == num_matches and hsp.query not in matches:
                    break
                matches.setdefault(hsp.query, []).append(hsp)
            if len(matches) > 0:
                seqs.append((hit, matches))
        return seqs
    def get_all_hits(self, max_e_value=None, min_identity=None, min_positive=None, min_query_coverage=None, min_hit_coverage=None, max_gap_coverage=None, min_hsp_length=None, min_bitscore=None, min_score=None, sort=None):
        """Returns a list of all query BlastSequence objects that have at least 1 HSP that satisfies all of the given filter criteria. The 'max_e_value', 'min_hsp_length', 'min_bitscore', and 'min_score' criteria all take float values; 'min_identity', 'min_positive', and 'max_gap_coverage' all take percent values (floats between 0 and 100) and are calculated over the length of the HSP; 'min_query_coverage' and 'min_hit_coverage' take percent values (floats between 0 and 100) indicating the minimum portion of the total query or hit length that must be aligned in an HSP. If 'sort' is not specified the returned list will be in the same order as the results file, otherwise it must be a single string from: "max_e_value", "min_identity", "min_positive", "min_query_coverage", "min_hit_coverage", "max_gap_coverage", "min_hsp_length", "min_bitscore", "min_score"."""
        hits = []
        for hit in self.hits.values():
            if hit.filter_hsps(max_e_value, min_identity, min_positive, min_query_coverage, min_hit_coverage, max_gap_coverage, min_hsp_length, min_bitscore, min_score):
                hits.append(hit)
        if sort == None:
            return hits
        return self.sort_sequences(hits, sort)
    def get_all_queries(self, max_e_value=None, min_identity=None, min_positive=None, min_query_coverage=None, min_hit_coverage=None, max_gap_coverage=None, min_hsp_length=None, min_bitscore=None, min_score=None, sort=None):
        queries = []
        for query in self.queries.values():
            if query.filter_hsps(max_e_value, min_identity, min_positive, min_query_coverage, min_hit_coverage, max_gap_coverage, min_hsp_length, min_bitscore, min_score):
                queries.append(query)
        if sort == None:
            return queries
        return self.sort_sequences(queries, sort)
    def sort_sequences(self, sequences, sort_key, filter_unaligned=False, sort_direction="best", hsp_value_type="best"):
        # sort_direction can be: "best", "worst", "ascending", or "descending". If "best", returns in best-to-worst order (lowest to highest e-value, highest to lowest identity, etc), reversed gives worst-to-best. hsp_value_type can be "best", "worst", "highest", "lowest", or "averge". If 'filter_unaligned' is True sequences that have no HSPs will be removed; if 'filter_unaligned' is False and 'sort_direction' is "worst" these will be at the front of the returned list, otherwise they will be at the end.
        seqs = [seq for seq in sequences if len(seq.hsps) > 0]
        no_hsps = [seq for seq in sequences if len(seq.hsps) == 0]
        if sort_direction == "best":
            if sort_key in self.best_is_low:
                seqs.sort(key=lambda seq: seq.get_hsp_value(sort_key, hsp_value_type))
            else:
                seqs.sort(key=lambda seq: seq.get_hsp_value(sort_key, hsp_value_type), reverse=True)
        elif sort_direction == "worst":
            if sort_key in self.best_is_low:
                seqs.sort(key=lambda seq: seq.get_hsp_value(sort_key, hsp_value_type), reverse=True)
            else:
                seqs.sort(key=lambda seq: seq.get_hsp_value(sort_key, hsp_value_type))
        elif sort_direction == "ascending":
            seqs.sort(key=lambda seq: seq.get_hsp_value(sort_key, hsp_value_type))
        elif sort_direction == "descending":
            seqs.sort(key=lambda seq: seq.get_hsp_value(sort_key, hsp_value_type), reverse=True)
        if filter_unaligned == True:
            return seqs
        if sort_direction == "worst":
            return no_hsps + seqs
        return seqs + no_hsps


    # # #  Parsing methods  # # #
    def parse_xml_data(self, xml_data):
        self.root = ET.fromstring(xml_data)
        if self.root.tag != 'BlastOutput':
            print('Error: unexpected file format. The root tag is "{}"'.format(self.root.tag))
            exit()
        self.blast_program = self.root.find('BlastOutput_version').text
        blast_type = self.root.find('BlastOutput_program').text
        self.parse_run_parameters()
        iterations = self.root.find('BlastOutput_iterations')
        if iterations == None:
            print('Error: unexpected file format. No "BlastOutput_iterations" tag found.')
            exit()
        self.queries = OrderedDict()
        self.hits = OrderedDict()
        for query_ele in iterations.findall('Iteration'):
            self.parse_iteration(query_ele)
    def parse_run_parameters(self):
        self.run_parameters = {}
        run_params = self.root.find('BlastOutput_param').find('Parameters')
        for tag, key in (('Parameters_matrix', 'matrix'), ('Parameters_expect', 'expect'), ('Parameters_gap-open', 'gap_open'), ('Parameters_gap-extend', 'gap_extend'), ('Parameters_filter', 'filter')):
            sub = run_params.find(tag)
            self.run_parameters[key] = sub.text if sub != None else None
    def parse_iteration(self, query_ele):
        q_id = query_ele.find('Iteration_query-ID').text
        if q_id in self.queries:
            print('Warning: a query sequence has a repeated identifier "{}". It will be ignored.'.format(q_id))
            return
        q_def = query_ele.find('Iteration_query-def').text
        q_acc = None
        q_len = query_ele.find('Iteration_query-len').text
        query = BlastSequence(self, q_id, q_def, q_acc, q_len)
        self.queries[q_id] = query
        for hit_ele in query_ele.find('Iteration_hits').findall('Hit'):
            hit_id = hit_ele.find('Hit_id').text
            if hit_id in self.hits:
                hit = self.hits[hit_id]
            else:
                hit_def = hit_ele.find('Hit_def').text
                hit_acc = hit_ele.find('Hit_accession').text
                hit_len = hit_ele.find('Hit_len').text
                hit = BlastSequence(self, hit_id, hit_def, hit_acc, hit_len)
                self.hits[hit_id] = hit
            for hsp_ele in hit_ele.find('Hit_hsps').findall('Hsp'):
                h_bitscore = hsp_ele.find('Hsp_bit-score').text
                h_score = hsp_ele.find('Hsp_score').text
                h_e_val = hsp_ele.find('Hsp_evalue').text
                q_range = (hsp_ele.find('Hsp_query-from').text, hsp_ele.find('Hsp_query-to').text)
                h_range = (hsp_ele.find('Hsp_hit-from').text, hsp_ele.find('Hsp_hit-to').text)
                h_idents = hsp_ele.find('Hsp_identity').text
                h_pos = hsp_ele.find('Hsp_positive').text
                h_gaps = hsp_ele.find('Hsp_gaps').text
                h_len = hsp_ele.find('Hsp_align-len').text
                q_seq = hsp_ele.find('Hsp_qseq').text
                h_seq = hsp_ele.find('Hsp_hseq').text
                hsp = BlastHsp(query, hit, h_e_val, h_idents, h_pos, q_range, h_range, h_gaps, h_len, h_bitscore, h_score, q_seq, h_seq)
                query.hsps.append(hsp)
                hit.hsps.append(hsp)

    # # #  Misc methods  # # #


class BlastSequence(object):
    """An object representing either a query or hit sequence"""
    def __init__(self, blast_results, id, definition, accession, length):
        self.blast_results = blast_results
        self.id = id
        self.definition = definition
        self.accession = accession
        self.length = int(length)
        self.description_str = '{} {}'.format(self.id, self.definition)
        if self.accession != None:
            self.description_str += ' (Accession: {})'.format(self.accession)
        self.hsps = []
    def filter_hsps(self, max_e_value=None, min_identity=None, min_positive=None, min_query_coverage=None, min_hit_coverage=None, max_gap_coverage=None, min_hsp_length=None, min_bitscore=None, min_score=None):
        filtered = []
        for hsp in self.hsps:
            if (max_e_value == None or hsp.e_value <= max_e_value) \
            and (min_identity == None or hsp.percent_identity >= min_identity) \
            and (min_positive == None or hsp.percent_positive >= min_positive) \
            and (min_query_coverage == None or hsp.query_coverage >= min_query_coverage) \
            and (min_hit_coverage == None or hsp.hit_coverage >= min_hit_coverage) \
            and (max_gap_coverage == None or hsp.gap_coverage <= max_gap_coverage) \
            and (min_hsp_length == None or hsp.hsp_length >= min_hsp_length) \
            and (min_bitscore == None or hsp.bit_score >= min_bitscore) \
            and (min_score == None or hsp.score >= min_score):
                filtered.append(hsp)
        return filtered
    def get_hsp_value(self, hsp_key, hsp_value_type="best"):
        """Returns a single value from the sequence's HSPs. The value retrieved is specified by the string 'hsp_key', which must be one of: "e_value", "percent_identity", "percent_positive", "query_coverage", "hit_coverage", "gap_coverage", "hsp_length", "bit_score", or "score". The string 'hsp_value_type' specifies which HSP value is returned and must be one of: "best" (the lowest e_value or gap_coverage, highest for the other keys), "worst", "highest", "lowest", or "average". If the sequence is associated with no HSPs, the worst possible value is returned (infinity for e_value and gap_coverage, zero for the other keys). This method is also used to sort sequences."""
        hsp_key = hsp_key.lower()
        if hsp_key not in self.blast_results.hsp_keys:
            print('Error: invalid key passed to BlastSequence.get_hsp_value() "{}".'.format(hsp_key))
            exit()
        hsp_vals = [getattr(hsp, hsp_key) for hsp in self.hsps]
        if not hsp_vals:
            if hsp_key in self.blast_results.best_is_low:
                return float('infinity')
            else:
                return 0.0
        if hsp_value_type == "best":
            if hsp_key in self.blast_results.best_is_low:
                return min(hsp_vals)
            else:
                return max(hsp_vals)
        elif hsp_value_type == "worst":
            if hsp_key in self.blast_results.best_is_low:
                return max(hsp_vals)
            else:
                return min(hsp_vals)
        elif hsp_value_type == "highest":
            return max(hsp_vals)
        elif hsp_value_type == "lowest":
            return min(hsp_vals)
        elif hsp_value_type == "average":
            return sum(hsp_vals) / float(len(hsp_vals))

class BlastHsp(object):
    """Object representing a high-scoring segment pair, which is an alignment between the query and hit."""
    def __init__(self, query, hit, e_value, identities, positives, query_range, hit_range, gaps, hsp_length, bit_score, score, query_sequence, hit_sequence):
        # #  BlastSequence references  # #
        self.query = query
        self.hit = hit
        # #  Raw values  # #
        self.e_value = float(e_value)
        self.identities = int(identities)
        self.positives = int(positives)
        self.query_range = (int(query_range[0]), int(query_range[1]))
        self.hit_range = (int(hit_range[0]), int(hit_range[1]))
        self.gaps = int(gaps)
        self.hsp_length = int(hsp_length)
        self.bit_score = float(bit_score)
        self.score = float(score)
        self.query_sequence = query_sequence
        self.hit_sequence = hit_sequence
        self.validate_values()
        # #  Calculated values  # #
        self.percent_identity = float(identities) / self.hsp_length * 100.0
        self.percent_positive = float(positives) / self.hsp_length * 100.0
        self.query_coverage = float(self.query_range[1] - self.query_range[0] + 1) / query.length * 100.0
        self.hit_coverage = float(self.hit_range[1] - self.hit_range[0] + 1) / hit.length * 100.0
        self.gap_coverage = float(gaps) / self.hsp_length * 100.0
    def validate_values(self):
        # ensure values are sane. If not, throw error.
        pass


# # #  Argument parser and option validation  # # #
def setup_parser():
    prog_descrip = 'A script to parse a BLAST results file in XML format, and perform some common functions. Note: the COMMAND argument must appear after all non-command options.'
    prog_epilog = '' # probably a couple example usages
    parser = argparse.ArgumentParser(description=prog_descrip, epilog=prog_epilog, add_help=False)
    add_arguments(parser)
    return parser
def add_arguments(parser):
    # #  Positional arguments  # #
    parser.add_argument("xml_file", metavar="XML_FILE", help="The BLAST results file in XML format")
    command_parser = parser.add_subparsers(dest="command", metavar="COMMAND", help="Must be one of:")
    command_parser.add_parser("summary", help="Get a summary description (default)")
    command_parser.add_parser("getids", help="Get all specified sequence IDs")
    command_parser.add_parser("getaccs", help="Get all specified sequence accessions")
    command_parser.add_parser("getdefs", help="Get all specified sequence definitions (note that queries usually have none)")
    command_parser.add_parser("getseqs", help="Get the sequence of all specified HSPs (only the aligned region, not the whole sequence)")
    # #  Sequence type arguments  # #
    seq_type_group = parser.add_argument_group(title="get results from", description="only one option may be specified")
    seq_type_ex_group = seq_type_group.add_mutually_exclusive_group()
    seq_type_ex_group.add_argument("-h", "--hits", action="store_true", help="Return information from the hits sequences (default)")
    seq_type_ex_group.add_argument("-q", "--queries", action="store_true", help="Return information from the query sequences")
    # #  Optional arguments  # #
    parser.add_argument("--help", action='help', help='show this help message and exit')
    parser.add_argument("-o", "--outfile", help="The name of a file to store the program output. If not given, the output will be printed")
    parser.add_argument("-r", "--range", action="store_true", help="Include the range of the match (only works with getids, getaccs, and getdefs)")
    parser.add_argument("-n", "--number_matches", type=int, help="Return information from this many sequences")
    parser.add_argument("-s", "--sort", choices=["e_value", "percent_identity", "percent_positive", "query_coverage", "hit_coverage", "gap_coverage", "hsp_length", "bit_score", "score"], nargs=1, metavar="SORT_KEY", help="Sort returned sequences; must be one of: %(choices)s")
    # #  Filtering arguments  # #
    filter_group = parser.add_argument_group(title="sequence filtering", description="these options allow the results to be filtered using one or more metrics")
    filter_group.add_argument("-e", "--e_value", type=float, help="Exclude matches with an e-value above '%(dest)s'")
    filter_group.add_argument("-i", "--identity", type=float, help="Exclude matches with a percent identity below '%(dest)s'")
    filter_group.add_argument("-p", "--positive", type=float, help="Exclude matches with a percent positive below '%(dest)s'")
    filter_group.add_argument("-qc", "--query_coverage", type=float, help="Exclude matches that cover less than '%(dest)s' percent of a query")
    filter_group.add_argument("-hc", "--hit_coverage", type=float, help="Exclude matches that cover less than '%(dest)s' percent of a hit")
    filter_group.add_argument("-gc", "--gap_coverage", type=float, help="Exclude matches where gaps cover more than '%(dest)s' percent of an HSP")
    filter_group.add_argument("-hl", "--hsp_length", type=int, help="Exclude matches with HSP alignments shorter than '%(dest)s'")
    filter_group.add_argument("-bs", "--bit_score", type=float, help="Exclude matches with an alignment bit score less than '%(dest)s'")
    filter_group.add_argument("-as", "--alignment_score", type=float, help="Exclude matches with an alignment score less than '%(dest)s'; it's generally better to specify the bit_score instead of this score")
def get_and_validate_arguments(parser):
    """Also modifies args.xml_file, args.hits, and args.sort"""
    args = parser.parse_args()
    # #  Positional arguments  # #
    args.xml_file = os.path.realpath(args.xml_file)
    if not os.path.isfile(args.xml_file):
        parser.error('unable to locate XML file at {}'.format(args.xml_file))
    # #  Sequence type arguments  # #
    if args.hits == args.queries == False:
        args.hits = True
    # #  Optional arguments  # #
    if args.range == True and args.command not in ("getids", "getaccs", "getdefs"):
        parser.error('-r/--range can only be used if the command is one of: getids, getaccs, getdefs')
    if args.number_matches != None and args.number_matches <= 0:
        parser.error('-n/--number must be an integer greater than 0')
    if args.sort != None:
        args.sort = args.sort[0]
    # #  Filtering arguments  # #
    if args.e_value != None and args.e_value < 0:
        parser.error('-e/--e_value must be greater than or equal to 0')
    if args.identity != None and not 0 < args.identity <= 100:
        parser.error('-i/--identity must be greater than 0 but less than or equal to 100')
    if args.positive != None and not 0 < args.positive <= 100:
        parser.error('-p/--positive must be greater than 0 but less than or equal to 100')
    if args.query_coverage != None and not 0 < args.query_coverage <= 100:
        parser.error('-qc/--query_coverage must be greater than 0 but less than or equal to 100')
    if args.hit_coverage != None and not 0 < args.hit_coverage <= 100:
        parser.error('-hc/--hit_coverage must be greater than 0 but less than or equal to 100')
    if args.gap_coverage != None and not 0 <= args.gap_coverage < 100:
        parser.error('-gc/--gap_coverage must be greater than or equal to 0 but less than 100')
    if args.hsp_length != None and args.hsp_length <= 0:
        parser.error('-hl/--hsp_length must be an integer greater than 0')
    if args.bit_score != None and args.bit_score <= 0:
        parser.error('-bs/--bit_score must be greater than 0')
    if args.alignment_score != None and args.alignment_score <= 0:
        parser.error('-as/--alignment_score must be greater than 0')
    return args


if __name__ == '__main__':
    try:
        import argparse
    except ImportError:
        print('Error: Python is unable to import the module "argparse", so this module cannot be run from the command line. This is most likely because you have an old version of Python.')
        exit()
    parser = setup_parser()
    args = get_and_validate_arguments(parser)

    buff = []
    results = results_from_file(args.xml_file)
    if args.command == "summary":
        max_line_width = 80 # Used for formatting output
        if args.hits == True:
            seq_info = results.get_hits_per_query(max_e_value=args.e_value, min_identity=args.identity, min_positive=args.positive, min_query_coverage=args.query_coverage, min_hit_coverage=args.hit_coverage, max_gap_coverage=args.gap_coverage, min_hsp_length=args.hsp_length, min_bitscore=args.bit_score, min_score=args.alignment_score, num_matches=args.number_matches, sort=args.sort)
        else:
            seq_info = results.get_queries_per_hit(max_e_value=args.e_value, min_identity=args.identity, min_positive=args.positive, min_query_coverage=args.query_coverage, min_hit_coverage=args.hit_coverage, max_gap_coverage=args.gap_coverage, min_hsp_length=args.hsp_length, min_bitscore=args.bit_score, min_score=args.alignment_score, num_matches=args.number_matches, sort=args.sort)
        for main_seq, matches in seq_info:
            buff.append(main_seq.description_str)
            buff.append("=" * min(len(main_seq.description_str), max_line_width))
            for match_seq, hsps in matches.items():
                buff.append(match_seq.description_str)
                buff.append("-" * min(len(match_seq.description_str), max_line_width))
                for hsp in hsps:
                    buff.append("- E-value {0.e_value:.2g}, Bit-score {0.bit_score:.1f} | Query {0.query_range[0]} - {0.query_range[1]} ({0.query_coverage:.1f}%), Hit {0.hit_range[0]} - {0.hit_range[1]} ({0.hit_coverage:.1f}%) | Identities {0.identities} ({0.percent_identity:.1f}%), Positives {0.positives} ({0.percent_positive:.1f}%) | Gaps {0.gaps} ({0.gap_coverage:.2f}%), HSP length {0.hsp_length}".format(hsp))
                buff[-1] += '\n'
            buff[-1] += '\n'
    elif args.command in ("getids", "getaccs", "getdefs", "getseqs"):
        if args.hits == True:
            seqs = results.get_all_hits(max_e_value=args.e_value, min_identity=args.identity, min_positive=args.positive, min_query_coverage=args.query_coverage, min_hit_coverage=args.hit_coverage, max_gap_coverage=args.gap_coverage, min_hsp_length=args.hsp_length, min_bitscore=args.bit_score, min_score=args.alignment_score, sort=args.sort)
        else:
            seqs = results.get_all_queries(max_e_value=args.e_value, min_identity=args.identity, min_positive=args.positive, min_query_coverage=args.query_coverage, min_hit_coverage=args.hit_coverage, max_gap_coverage=args.gap_coverage, min_hsp_length=args.hsp_length, min_bitscore=args.bit_score, min_score=args.alignment_score, sort=args.sort)
        if args.command == "getseqs":
            # for a seq, need to find the best (longest? bitscore?) hsp, and just take the aligned seq from that. getseqs should also have a specific option to remove gaps or not.
            # Note: stopcodons (*) are replaced with an X
            for seq in seqs:
                best_len, best_seq = 0, None
                for hsp in seq.hsps:
                    if hsp.hsp_length > best_len:
                        best_len = hsp.hsp_length
                        best_seq = hsp.hit_sequence if args.hits else hsp.query_sequence
                if best_seq:
                    buff.append('>{}'.format(seq.description_str))
                    buff.append('{}\n'.format(best_seq.replace('-','').replace('*','X')))
            print('Got {} sequences.'.format(len(seqs)))
        else:
            if args.command == "getids":
                seq_attr = 'id'
                command_type = 'IDs'
            elif args.command == "getaccs":
                seq_attr = 'accession'
                command_type = 'accessions'
            elif args.command == "getdefs":
                seq_attr = 'definition'
                command_type = 'definitions'
            range_attr = 'hit_range' if args.hits == True else 'query_range'
            for seq in seqs:
                seq_desc = getattr(seq, seq_attr, None)
                if seq_desc == None:
                    continue
                buff.append(seq_desc)
                if args.range == True:
                    ranges = sorted(set(getattr(hsp, range_attr) for hsp in seq.hsps))
                    buff.append('\n'.join('  {0[0]} - {0[1]}'.format(rng) for rng in ranges))
            print('Got {} from {} sequences.'.format(command_type, len(seqs)))

    buff_str = '\n'.join(buff) + '\n' # Ends in a newline.
    if args.outfile != None:
        out_path = os.path.realpath(args.outfile)
        with open(out_path, 'w') as f:
            f.write(buff_str)
        print('Output saved to {}'.format(out_path))
    else:
        print(buff_str)
