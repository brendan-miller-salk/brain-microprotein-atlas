import argparse
import file_reader
import statistics


def parse_args():
    parser = argparse.ArgumentParser(
        description=('Normalize the raw counts in the abundance file'))
    parser.add_argument('--abundance-esp',
                        help='the *_abundance.esp file output by ESPRESSO',
                        required=True)
    parser.add_argument('--output-path',
                        help='base path for saving output files (unfiltered, filtered, and ids.txt)',
                        required=True)

    return parser.parse_args()


def get_totals_by_sample(in_esp_path):
    with open(in_esp_path, 'rt') as in_esp:
        for line_i, line in enumerate(in_esp):
            if line_i == 0:
                initial_headers, sample_headers = (
                    file_reader.read_abundance_esp_header_line(line))
                sample_totals = [0] * len(sample_headers)
                continue

            values = file_reader.read_abundance_esp_line(line, sample_headers)
            for i, sample in enumerate(sample_headers):
                sample_totals[i] += values[sample]

    return sample_totals


def write_columns(out_f, columns):
    str_columns = list()
    for col in columns:
        if isinstance(col, float):
            str_val = '{:.2f}'.format(col)
        else:
            str_val = str(col)

        str_columns.append(str_val)

    out_f.write('{}\n'.format('\t'.join(str_columns)))


def write_normalized_esp(totals_by_sample, in_esp_path, output_base_path):
    filtered_output_path = f"{output_base_path}_medianCPM05.txt"
    ids_output_path = f"{output_base_path}__medianCPM05_ids.txt"
    unfiltered_output_path = f"{output_base_path}_allCPM.txt"

    transcript_ids = []

    with open(in_esp_path, 'rt') as in_esp:
        with open(unfiltered_output_path, 'wt') as unfiltered_out, \
                open(filtered_output_path, 'wt') as filtered_out, \
                open(ids_output_path, 'wt') as ids_out:

            for line_i, line in enumerate(in_esp):
                if line_i == 0:
                    initial_headers, sample_headers = (
                        file_reader.read_abundance_esp_header_line(line))
                    write_columns(unfiltered_out, initial_headers + sample_headers)
                    write_columns(filtered_out, initial_headers + sample_headers)
                    continue

                values = file_reader.read_abundance_esp_line(
                    line, sample_headers)
                cpms = []
                for sample_i, sample in enumerate(sample_headers):
                    sample_val = values[sample]
                    total = totals_by_sample[sample_i]
                    cpms.append((sample_val * 1e6) / total)

                # Calculate median CPM
                median_cpm = statistics.median(cpms)

                # Write unfiltered data
                initial_columns = [values[h] for h in initial_headers]
                write_columns(unfiltered_out, initial_columns + cpms)

                # Write filtered data and collect transcript IDs if median CPM >= 0.5
                if median_cpm >= 0.5:
                    write_columns(filtered_out, initial_columns + cpms)
                    transcript_ids.append(initial_columns[0])  # Assuming the first column is the transcript ID

            # Save transcript IDs to ids.txt
            for transcript_id in transcript_ids:
                ids_out.write(f"{transcript_id}\n")


def main():
    args = parse_args()
    totals_by_sample = get_totals_by_sample(args.abundance_esp)
    write_normalized_esp(totals_by_sample, args.abundance_esp, args.output_path)


if __name__ == '__main__':
    main()