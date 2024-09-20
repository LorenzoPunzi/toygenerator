import re
import os
import matplotlib
matplotlib.use('Agg')  # Use the 'Agg' backend to avoid opening windows
import matplotlib.pyplot as plt

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
            data_match = re.match(r"\s*([-\d.]+)\s+([-\d.]+)\s+([\d.]+)\s+([\d.]+)", line)
            if data_match and current_histogram:
                binlow = float(data_match.group(1))
                binhigh = float(data_match.group(2))
                content = float(data_match.group(3))
                error = float(data_match.group(4))
                histograms[current_histogram].append((binlow, binhigh, content, error))

    return histograms

def plot_histogram(hist_name, bins, output_dir):
    # Prepare the data for plotting
    bin_edges = [binlow for binlow, binhigh, _, _ in bins] + [bins[-1][1]]  # Add the last bin's upper edge
    contents = [content for _, _, content, _ in bins]
    errors = [error for _, _, _, error in bins]

    # Create the plot
    plt.figure(figsize=(8, 6))
    plt.bar(bin_edges[:-1], contents, width=[binhigh - binlow for binlow, binhigh, _, _ in bins], 
            align='edge', yerr=errors, capsize=5, alpha=0.7, color='b', label='Content')

    # Add labels and title
    plt.xlabel('Bins')
    plt.ylabel('Content')
    plt.title(f'Histogram: {hist_name}')

    # Save the plot as a PNG file in the output directory
    plt.tight_layout()
    output_path = os.path.join(output_dir, f'{hist_name}.png')
    plt.savefig(output_path)
    plt.close()

def save_histograms(histograms, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    for hist_name, bins in histograms.items():
        plot_histogram(hist_name, bins, output_dir)

if __name__ == "__main__":
    # Replace 'hist.out' with the path to your histogram file
    file_path = 'hist.out'
    output_dir = 'histogram_plots'  # Directory to save the plots
    
    histograms = read_histogram(file_path)
    save_histograms(histograms, output_dir)

    print(f"Histograms have been saved to the directory: {output_dir}")
