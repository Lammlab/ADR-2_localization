"""Fastq Collapse

    Usage:
        fastq_collapse <fastq_filename> <new_filename> [ <prefix>] [--count]
        fastq_collapse param
        fastq_collapse example
        fastq_collapse -h | --help

    Options:
        -h --help   Show this screen
        --count     include a count of each gene in the collapsed file
"""

########################################################################################################################
# Refactor By:          Amit Mark Finkelstein (github: amitfinkel)
# Main goal:            The script takes a path to a fastq input file and a name or path (chosen by the user)
#                       for an output file, and returns a collapsed output file.
########################################################################################################################

import os
import sys
from docopt import docopt
from tqdm import tqdm
from Utility.Fastq_class import Fastq
from Utility.generators_utilities import BAR_DEFAULT_VIEW, generatesKLines
from Utility.generators_utilities import num_of_lines_in_file as calc_num_lines_in_file
DEFAULT_PREFIX = 3


def maximum_score(curr_score, dict_score):
    """
    :param curr_score: First score string to compare
    :param dict_score: Second score string to compare
    :return: A new string of score, where each index is the maximum score
            between the 2 score strings
    """
    new_score = ""
    for i in range(len(curr_score)):  # For each char in score
        curr_char = curr_score[i]
        dict_char = dict_score[i]
        new_score += curr_char if curr_char > dict_char else dict_char  # Put maximum

    return new_score


def collapse_fastq_to_dict(generator_file):
    """
    :param generator_file: A generator of 4 lines of the file of 'FastQ' elements'.
    :return: A dictionary of the format: { 'string' : ('FastQ' object, 'int') }. With 'string' being
             the base-pair seq, 'object' being the best 'FastQ' seq found,
             and 'int' being that particular base-pair seq counter
    """
    seq_dict = {}

    for lines in generator_file:
        last_line = lines[3] if lines[3][-1:] != "\n" else lines[3][:-1]
        seq = Fastq([lines[0][:-1], lines[1][:-1], lines[2][:-1], last_line])
        try:
            seq_dict[seq.strings[1]][1] += 1
            # if qualityIsBetter(seq.strings[3], seq_dict[seq.strings[1]][0].strings[3]):
            #    seq_dict[seq.strings[1]][0] = seq
            curr_score = seq.strings[3]
            dict_score = seq_dict[seq.strings[1]][0].strings[3]
            seq_dict[seq.strings[1]][0].strings[3] = maximum_score(curr_score, dict_score)
        except Exception as exception:  # This is the first occurrence of this sequence
            if type(exception) == KeyError:
                seq_dict[seq.strings[1]] = [seq, 1]
            else:
                print(exception)
                sys.exit(-1)

    return seq_dict


def generate_fastq_file_from_dict(fastq_dict, filename, generating_fastq_progress_bar, write_count_flag):
    """
    :param fastq_dict: The dictionary described above (See 'collapseFastqSeqListToDict.__doc__')
    :param filename: The requested filename (with relative path) of the new collapsed fastq file
    :param generating_fastq_progress_bar: A progress bar that will update while the new fastq file is being generated.
    :param write_count_flag: a flag indicating whether to write sequence counts to the output file
    :return: The same filename entered at the input OR -1 if a file with the same name exits
             in the designated directory.
    """
    if filename.split("/")[-1] in "".join(filename.split("/")[:-1]):
        print("A file with that name already exits in the designated directory!")
        return -1
    new_fastq_file = open(filename, "a+")  # Read + Append, 'line buffering'

    for fastq_obj, count in fastq_dict.values():
        if write_count_flag:
            fastq_obj.strings[0] = fastq_obj.strings[0] + "%s:%d" % ('_count', count)
        new_fastq_file.write("\n".join(fastq_obj.strings) + "\n")
        generating_fastq_progress_bar.update(1)

    os.chmod(filename, 0o777)
    new_fastq_file.close()

    return filename


def split_file_to_sub_files_by_prefix(generator_file, prefix, fastq_filename, generated_filename):
    """
    :param generator_file: A generator of 4 lines of the file of 'FastQ' elements
    :param prefix: A length for the prefix of nucleotides for the separation of the input fastq file into subfiles.
    :param fastq_filename: The requested filename (with relative path) of the fastq file that will be collapsed.
    :param generated_filename: The fastq filename for the generator and for tqdm usage.
    :return: a list of files names
    """
    list_of_file_names = []
    splitting_progress_bar = tqdm(total=calc_num_lines_in_file(generated_filename),
                                  desc="Splitting files into subfiles by prefixes ", bar_format=BAR_DEFAULT_VIEW)
    lines = generator_file.__next__()
    splitting_progress_bar.update(4)
    if not lines[0]:
        open(fastq_filename, 'w').close()
    # created a small (10000 lines from file for example) dictionary and push
    # to the fitted files by the prefix
    writing_progress_bar = tqdm(total=calc_num_lines_in_file(generated_filename),
                                desc="Writing data into temp files ", bar_format=BAR_DEFAULT_VIEW)
    writing_progress_bar.update(4)
    while lines[0]:
        try:
            dict_of_file_names = {}  # also frees the space after each time
            if not lines[0]:
                break
            first_letter_of_seq = lines[1][:prefix]
            file_name_for_seq = fastq_filename + "_" + first_letter_of_seq + ".temp"
            # so we could run more than one collapse at a time
            if file_name_for_seq not in list_of_file_names:
                list_of_file_names.append(file_name_for_seq)
            if file_name_for_seq not in dict_of_file_names:
                dict_of_file_names[file_name_for_seq] = [lines]
            else:
                dict_of_file_names[file_name_for_seq].append(lines)
            lines = generator_file.__next__()
            splitting_progress_bar.update(4)

        except StopIteration:
            for name in dict_of_file_names.keys():
                with open(name, 'a+') as file:
                    full_data_for_file = ""
                    for lines_chunk in dict_of_file_names[name]:
                        original_chunk = "".join(lines_chunk)
                        full_data_for_file = full_data_for_file + original_chunk
                    file.write(full_data_for_file)
            break

        for name in dict_of_file_names.keys():
            with open(name, 'a+') as file:
                full_data_for_file = ""
                for lines_chunk in dict_of_file_names[name]:
                    original_chunk = "".join(lines_chunk)
                    full_data_for_file = full_data_for_file + original_chunk
                file.write(full_data_for_file)
            writing_progress_bar.update(4)

    return list_of_file_names


def fastq_collapse(fastq_filename, new_filename, prefix=3, write_count_flag=True):
    """
    :param fastq_filename: The requested filename (with relative path) of the fastq file that will be collapsed.
    :param new_filename: The requested filename (with relative path) of the fastq file that will be created.
    :param prefix: A length for the prefix of nucleotides for the separation of the input fastq file into subfiles.
    :param write_count_flag: a flag indicating whether to write sequence counts to the output file
    :return: a list of files names
        """
    with open(fastq_filename, "r") as fastq_file:  # Read + Write, 'line buffering'
        # creates a generator for the file
        generator_file = generatesKLines(fastq_file, num_lines=4)
        # splitting the fastqs to files by their seq prefix
        list_of_files = split_file_to_sub_files_by_prefix(generator_file, prefix, new_filename, fastq_filename)
        # collapsing each file
        generating_fastq_progress_bar = tqdm(desc="Generating Fastq file from newly created dictionary ",
                                             bar_format=BAR_DEFAULT_VIEW)
        for filename in tqdm(list_of_files,
                             desc="Collapsing into relevant files ", bar_format=BAR_DEFAULT_VIEW):
            # creating generator for file
            with open(filename, "r") as file:
                generator_file = generatesKLines(file, num_lines=4)
                # collapse each file
                collapsed_fastq_dict = collapse_fastq_to_dict(generator_file)  # "collapsing" it (to a dict)
                generating_fastq_progress_bar.reset(total=len(collapsed_fastq_dict))
                # sending all the collapsing to the same file
                new_filename = generate_fastq_file_from_dict(collapsed_fastq_dict, new_filename,
                                                             generating_fastq_progress_bar, write_count_flag)
                # generating the new "collapsed" file
            os.remove(filename)


def param_description():
    """
    Description for the script
    """
    print("The parameters are\n" +
          "fastq_filename: the name of the input fastq file\n"
          "new_filename: the name of the new file that would be created\n"
          "numLinesIter: size of chunks to read to memory \n"
          "prefix: size of prefix we split the files by while processing \n"
          "count: optional flag, indicating whether the collapsed file should include counts of each sequence"
          "Output: \"collapsed\" fastq format file")


def example_description():
    """
    Usage example for the script
    """
    print("Fastq file input:\n" +
          "@ABC\nGATC\n+\n!#!!\n" +
          "@BCD\nGATC\n+\n\"\"\"\"\n" +
          "@GHJ\nTCGA\n+\n####\n\n" +
          "Calling the script with count flag:\n" +
          "python fastq_collapse.py fastq_input fastq_output --count\n\n"
          "Output file:\n" +
          "@ABC_count:2\nGATC\n+\n!\"!!\n" +
          "@GHJ_count:1\nTCGA\n+\n####\n\n" +
          "Calling the script without count flag:\n" +
          "python fastq_collapse.py fastq_input fastq_output\n\n"
          "Output file:\n" +
          "@ABC\nGATC\n+\n!\"!!\n" +
          "@GHJ\nTCGA\n+\n####\n")


# ********************************************************

if __name__ == "__main__":

    """
    gets a fastq file , a newName for output file , a number of lines to iterates each time ,
     number of letters to split files by
    """

    arguments = docopt(__doc__)

    if arguments["param"]:
        param_description()
        sys.exit()

    if arguments["example"]:
        example_description()
        sys.exit()

    try:
        prefix = DEFAULT_PREFIX
        write_count_flag = True
        fastq_filename = arguments["<fastq_filename>"]
        new_filename = arguments["<new_filename>"]
        if not arguments["--count"]:
            write_count_flag = False
        if arguments["<prefix>"]:
            prefix = arguments["<prefix>"]
        fastq_collapse(fastq_filename, new_filename, prefix, write_count_flag)

    except Exception as exp:
        print(exp)
        sys.exit(2)

    with open("finished_actions.log", "a+") as logfile:
        logfile.write(arguments["<new_filename>"] + " finished")
# ********************************************************
