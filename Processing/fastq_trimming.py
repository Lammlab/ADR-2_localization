"""Fastq Trimming

    Usage:
        fastq_trimming <fastq_filename> <new_filename> <range_start> <range_end>
        fastq_trimming param
        fastq_trimming example
        fastq_trimming -h | --help

    Options:
        -h --help   Show this screen
"""

########################################################################################################################
# Refactor By:          Amit Mark Finkelstein (github: amitfinkel)
# Main goal:            The script takes a path to a fastq input file and a name or path (chosen by the user)
#                       for an output file, and trims the fastq sequences from the positions given by the user
#                       (please note: the positions should be zero-based).
#                       For example: trimming of sequence 'AAGGCT' with range_start=2 and range_end=5 will give
#                       the following sequence: 'GGC'.
########################################################################################################################

import sys
from docopt import docopt
from Utility.Fastq_class import Fastq
from Utility.generators_utilities import generatesKLines
from tqdm import tqdm
from Utility.generators_utilities import BAR_DEFAULT_VIEW
from Utility.generators_utilities import num_of_lines_in_file


def trimmByRange(fastq_filename, out_filename, start, end):
    fastq_file = open(fastq_filename, "r")
    lines = generatesKLines(fastq_file, num_lines=4)
    reading_progress_bar = tqdm(total=num_of_lines_in_file(fastq_filename), desc='Reading Input ',
                                bar_format=BAR_DEFAULT_VIEW)
    list_to_print = []
    last_line = ""
    for chunk in lines:
        if not chunk[0] or chunk is None:
            break
        last_line = chunk[3] if chunk[3][-1:] != "\n" else chunk[3][:-1]
        seq = Fastq([chunk[0][:-1], chunk[1][:-1], chunk[2][:-1], last_line])
        new_seq = seq.cut_seq(start, end)
        list_of_lines = new_seq.strings
        chunk_of_lines = ""
        for line in list_of_lines:
            if line[-1:] != "\n":
                chunk_of_lines = chunk_of_lines + line + "\n"
            else:
                chunk_of_lines = chunk_of_lines + line
        list_to_print.append(chunk_of_lines)
        reading_progress_bar.update(4)
    with open(out_filename, "w") as out_fp:
        for string in tqdm(list_to_print, desc='Writing To File ', bar_format=BAR_DEFAULT_VIEW):
            out_fp.write(string)
    fastq_file.close()


def fastq_trimming(fastq_filename, new_filename, range_start, range_end):
    if type(int(range_start)) != int or type(int(range_end)) != int:
        print("Incorrect entered start and end of the range - should be numbers")
        sys.exit(2)

    trimmByRange(fastq_filename, new_filename, int(range_start), int(range_end))


def param_description():
    """
    Description for the script
    """
    print("The parameters are\n" +
          "fastq_filename: the name of the input fastq file\n" +
          "new_filename: the name of the new file that would be created\n" +
          "range_start: start of trimming\n" +
          "end_start: end of trimming\n" +
          "Output: trimmed fastq format file")


def example_description():
    """
    Usage example for the script
    """
    print("Fastq file input:\n" +
          "@ABC\nGATCATC\n+\n!#!!$!@\n" +
          "@BCD\nGATCTTA\n+\n\"\"\"\"!!!\n" +
          "@GHJ\nTCGAGCAG\n+\n####!!!!\n\n" +
          "Calling the script:\n" +
          "python fastq_trimming fastq_input fastq_trimmed 2 5\n\n" +
          "Output file:\n" +
          "@ABC\nTCA\n+\n!!$\n" +
          "@BCD\nTCT\n+\n\"\"!\n" +
          "@GHJ\nGAG\n+\n##!")


if __name__ == "__main__":

    """
    gets a fastq file , a newName for output file , a start range, end range , num lines of file to iterate each time
    """

    arguments = docopt(__doc__)

    if arguments["param"]:
        param_description()
        sys.exit()

    if arguments["example"]:
        example_description()
        sys.exit()

    fastq_filename = arguments["<fastq_filename>"]
    new_filename = arguments["<new_filename>"]
    range_start = arguments["<range_start>"]
    range_end = arguments["<range_end>"]

    try:
        if fastq_filename and new_filename and range_start and range_end:
            fastq_trimming(fastq_filename, new_filename,
                           range_start, range_end)
    except Exception as exp:
        print(exp)
        sys.exit(2)
