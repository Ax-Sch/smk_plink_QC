import argparse
import re

def update_custom_log(input_log_file, custom_logs, output_file_prefix, custom_info=None):
"""Updates the custom log file with extracted information from the tool's log file and optional custom information."""
    with open(input_log_file, 'r') as infile, open(custom_logs, 'a') as outfile:
        outfile.write(f"Output File: {output_file_name}\n")
        outfile.write("-------------------------------\n")
        for line in infile:
            if re.search("Total genotyping rate", line) or \
            re.search("variants removed due to missing genotype", line) or \
            re.search("variants and people pass filters", line) or \
            re.search("cases and controls", line)
            outfile.write(line)

        if custom_info:
        outfile.write(f"Custom Information: {custom_info}\n")
        outfile.write("-------------------------------\n")


    if __name__ == "__main__":
        parser = argparse.ArgumentParser(description="Update custom log file with extracted information.")
        parser.add_argument("--input_log_file",required=True, help="Path to the produced log file")
        parser.add_argument("--custom_logs", required=True, help="Path to the custom log file")
        parser.add_argument("--output_file_prefix", required=True, help="prefix of the output file")
        parser.add_argument("--custom_info", help="Optional custom information", default=None)
                     
        # Parse arguments
        args = parser.parse_args()

        # Call the function with the parsed arguments
        update_custom_log(args.input_log_file, args.custom_logs, args.output_file_prefix, args.custom_info)
