import sys
import numpy as np
import tempfile
import os

def process_bedgraph(input_file, output_file, chunk_size=100000):
    total_sum = 0
    total_rows = 0

    # Temporary file for sorted data
    temp_sorted_file = tempfile.NamedTemporaryFile(delete=False, mode="w")

    # Step 1: Read, sort, and write sorted chunks to a temporary file
    data_chunks = []
    with open(input_file, "r") as infile:
        chunk = []
        for line in infile:
            fields = line.strip().split("\t")
            chr_, start, end, coverage = fields[0], int(fields[1]), int(fields[2]), float(fields[3])
            chunk.append((chr_, start, end, coverage))
            total_sum += coverage
            total_rows += 1

            if len(chunk) >= chunk_size:
                data_chunks.append(np.array(chunk, dtype=object))
                chunk = []

        if chunk:  # Process any remaining data
            data_chunks.append(np.array(chunk, dtype=object))

    # Sort and write sorted data
    sorted_data = np.concatenate(data_chunks)
    sorted_data = sorted_data[np.lexsort((sorted_data[:,1].astype(int), sorted_data[:,0], sorted_data[:,3].astype(float)))]
    
    with open(temp_sorted_file.name, "w") as outfile:
        for row in sorted_data:
            outfile.write("\t".join(map(str, row)) + "\n")

    # Step 2: Process the sorted file line-by-line for cumulative sum and slope
    cumsum = 0
    prev_cumsum = None

    with open(temp_sorted_file.name, "r") as infile, open(output_file, "w") as outfile:
        for rank, line in enumerate(infile, start=1):
            fields = line.strip().split("\t")
            chr_, start, end, coverage = fields[0], int(fields[1]), int(fields[2]), float(fields[3])

            cumsum += coverage
            norm_cumsum = (cumsum * total_rows) / total_sum if total_sum > 0 else 0

            # Compute slope (first entry has no previous cumsum)
            slope = "NA" if prev_cumsum is None else norm_cumsum - prev_cumsum

            # Write output
            outfile.write(f"{chr_}\t{start}\t{end}\t{coverage}\t{rank}\t{norm_cumsum:.12f}\t{slope}\n")

            # Track previous cumsum for slope calculation
            prev_cumsum = norm_cumsum

    # Remove temporary sorted file
    os.remove(temp_sorted_file.name)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 process_bedgraph.py <input_bedgraph> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    process_bedgraph(input_file, output_file)
