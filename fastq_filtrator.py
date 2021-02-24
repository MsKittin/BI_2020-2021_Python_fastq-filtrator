"""
Program for filtering .fastq files

Authors: Ekaterina Kirillova (parsing)
         Maksim Serdakov (filtering)

Supported arguments:
Optional:
         --min_length INT
         Minimum length of reads to pass the filtering.

         --keep_filtered
         If specified - failed-filtering reads will be written to a separate file.

         --gc_bounds INT (INT)
         Takes from 1 to 2 values:
         if 1 value is specified - this is the lower bound of GC-content in read to pass the filtering,
         if 2 values are specified - these are the lower and upper bounds respectively in read to pass the filtering.

         --output_base_name STR
         Common prefix for output file(s),
         if not specified - output file(s) will be named as source file without extension.

Required:
         .fastq file
         Positional argument, must be specified last.
"""

import os
import sys


class Config:
    def __init__(self, output_base_name: str) -> None:
        self.min_length = None
        self.keep_filtered = False
        self.output_base_name = output_base_name
        self.min_gc_bound = None
        self.max_gc_bound = None


def init_config(config: Config, arguments: list):

    arg_val_index = 0

    while arg_val_index <= len(args) - 1:

        if arguments[arg_val_index] == "--min_length":
            if arg_val_index + 1 < len(arguments):
                min_length_value = arguments[arg_val_index + 1]

                if not min_length_value.isdigit():
                    raise TypeError('--min_length value must be an integer')
                if int(min_length_value) <= 0:
                    raise ValueError('--min_length value must be > 0')

                config.min_length = int(min_length_value)

                arg_val_index += 2
            else:
                raise ValueError('Value for --min_length expected')

        elif arguments[arg_val_index] == "--keep_filtered":
            config.keep_filtered = True

            arg_val_index += 1

        elif arguments[arg_val_index] == "--gc_bounds":
            if arg_val_index + 1 < len(arguments):
                gc_bounds_value_min = arguments[arg_val_index + 1]
                max_bounds_flag = False

                if (arg_val_index + 2) < len(arguments) and arguments[arg_val_index + 2].isdigit():
                    gc_bounds_value_max = arguments[arg_val_index + 2]

                    if int(gc_bounds_value_max) < 0:
                        raise ValueError('Value for upper gc bound must be >= 0')
                    if int(gc_bounds_value_max) < int(gc_bounds_value_min):
                        raise ValueError('Lower gc bound must be less than upper bound')

                    config.max_gc_bound = int(gc_bounds_value_max)

                    max_bounds_flag = True

                if not gc_bounds_value_min.isdigit():
                    raise TypeError('Lower gc bound must be an integer')
                if int(gc_bounds_value_min) < 0:
                    raise ValueError('Value for lower gc bound must be >= 0')

                config.min_gc_bound = int(gc_bounds_value_min)

                if max_bounds_flag:
                    arg_val_index += 3
                else:
                    arg_val_index += 2

            else:
                raise ValueError('Value for --gc_bounds expected')

        elif arguments[arg_val_index] == "--output_base_name":
            if arg_val_index + 1 < len(arguments):
                config.output_base_name = arguments[arg_val_index + 1]
                arg_val_index += 2

            else:
                raise ValueError('Value for --output_base_name expected')

        else:
            raise NameError('Unexpected argument: ' + arguments[arg_val_index])


def GC_percent(cur_read):
    return (cur_read.upper().count('C') + cur_read.upper().count('G')) / len(cur_read) * 100


def GC_bounds(seq, lower_bound, upper_bound):
    cur_percent = GC_percent(seq)
    if (cur_percent >= lower_bound) and (cur_percent <= upper_bound):
        return True
    return False


def min_len_find(seq, minlen):
    if len(seq) >= minlen:
        return True
    return False


def features_check(user_minlen, user_lower_bound, user_upper_bound):
    if user_lower_bound is None:
        work_lower_bound = 0
        work_upper_bound = 100
    elif config.max_gc_bound is None:
        work_lower_bound = user_lower_bound
        work_upper_bound = 100
    else:
        work_lower_bound = user_lower_bound
        work_upper_bound = user_upper_bound
    if user_minlen is None:
        work_minlen = False
    else:
        work_minlen = user_minlen
    return [work_lower_bound, work_upper_bound, work_minlen]


args = sys.argv[1:]

if not os.path.exists(args[-1]):
    raise FileExistsError('File does not exist: ' + args[-1])
if not os.path.isfile(args[-1]):
    raise IsADirectoryError('File as positional argument expected')


fastq_file_path = args[-1]
args = args[:-1]

arguments = args
config = Config(fastq_file_path.rsplit('.', 1)[0])
init_config(config, arguments)


working_parameters = features_check(config.min_length, config.min_gc_bound, config.max_gc_bound)
working_lower_bound = working_parameters[0]
working_upper_bound = working_parameters[1]
working_minlen = working_parameters[2]

if config.keep_filtered:
    good_reads_file = config.output_base_name + '__passed.fastq'
    bad_reads_file = config.output_base_name + '__failed.fastq'

    with open(fastq_file_path, 'r') as fastq_file, open(good_reads_file, 'w') as gf, open(bad_reads_file, 'w') as bf:
        full_data = [line.rstrip() for line in fastq_file]
        for i in range(1, len(full_data), 4):
            if working_minlen:
                min_len_result = True
            else:
                min_len_result = min_len_find(full_data[i], working_minlen)
            if GC_bounds(full_data[i], working_lower_bound, working_upper_bound) and min_len_result:
                gf.write(full_data[i-1] + '\n')
                gf.write(full_data[i] + '\n')
                gf.write(full_data[i+1] + '\n')
                gf.write(full_data[i+2] + '\n')
            else:
                bf.write(full_data[i - 1] + '\n')
                bf.write(full_data[i] + '\n')
                bf.write(full_data[i + 1] + '\n')
                bf.write(full_data[i + 2] + '\n')

else:
    good_reads_file = config.output_base_name

    with open(fastq_file_path, 'r') as fastq_file, open(good_reads_file, 'w') as gf:
        full_data = [line.rstrip() for line in fastq_file]
        for i in range(1, len(full_data), 4):
            if working_minlen:
                min_len_result = True
            else:
                min_len_result = min_len_find(data[i], working_minlen)
            if GC_bounds(full_data[i], working_lower_bound, working_upper_bound) and min_len_result:
                gf.write(full_data[i - 1] + '\n')
                gf.write(full_data[i] + '\n')
                gf.write(full_data[i + 1] + '\n')
                gf.write(full_data[i + 2] + '\n')
