import re

def read_histogram(file_path):
    histograms = {}
    current_histogram = None

    with open(file_path, 'r') as file:
        for line in file:
            # Match histogram names, like 'h_pmod1', 'h_emu1', etc.
            match = re.match(r"^\s*(h_\w+)", line)
            if match:
                current_histogram = match.group(1)
                histograms[current_histogram] = []
                continue
            
            # Match histogram data rows with binlow, binhigh, content, and error
            data_match = re.match(r"\s*([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)", line)
            if data_match and current_histogram:
                binlow = float(data_match.group(1))
                binhigh = float(data_match.group(2))
                content = float(data_match.group(3))
                error = float(data_match.group(4))
                histograms[current_histogram].append((binlow, binhigh, content, error))

    return histograms


def print_vertical_histogram(hist_name, bins, max_height=20, col_width=5):
    print(f"\nHistogram: {hist_name}")
    
    max_content = max(content for _, _, content, _ in bins) if bins else 1
    
    # Scale contents to fit within max_height
    scaled_bins = [int((content / max_content) * max_height) for _, _, content, _ in bins]
    
    # Transpose the histogram, print each row from top to bottom
    for row in range(max_height, 0, -1):
        line = ''
        for bin_height in scaled_bins:
            if bin_height >= row:
                line += ' # '.center(col_width)
            else:
                line += '   '.center(col_width)
        print(line)

    # Print the bin separators and labels
    print(" " + ("-" * col_width * len(bins)))
    
    # Print bin ranges below the bars, adjusting width for neatness
    for binlow, binhigh, _, _ in bins:
        bin_range = f"[{binlow:.2f},{binhigh:.2f}]"
        print(f"{bin_range:^{col_width}}", end=" ")
    print()


def print_histograms(histograms, col_width=7):
    for hist_name, bins in histograms.items():
        print_vertical_histogram(hist_name, bins, col_width=col_width)


if __name__ == "__main__":
    # Replace 'histogram_file.txt' with the path to your histogram file
    file_path = 'hist.out'
    histograms = read_histogram(file_path)
    print_histograms(histograms)
