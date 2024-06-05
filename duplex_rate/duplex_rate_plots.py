import matplotlib.pyplot as plt
import argparse
import os

# Function to extract data from the log content
def extract_data(log_content):
    lines = [line.split("[info]")[1]for line in log_content.splitlines()]
    #print(lines)
    simplex_reads_basecalled = 0
    simplex_reads_filtered = 0
    duplex_reads_basecalled = 0
    duplex_rate = 0.0

    for line in lines:
        if "Simplex reads basecalled" in line:
            simplex_reads_basecalled = int(line.split(":")[1].strip())
            #print('simplex_reads_basecalled:', simplex_reads_basecalled)
        elif "Simplex reads filtered" in line:
            simplex_reads_filtered = int(line.split(":")[1].strip())
            #print('simplex_reads_filtered:', simplex_reads_filtered)
        elif "Duplex reads basecalled" in line:
            duplex_reads_basecalled = int(line.split(":")[1].strip())
            #print('duplex_reads_basecalled:', duplex_reads_basecalled)
        elif "Duplex rate" in line:
            duplex_rate = float(line.split(":")[1].strip().replace('%', ''))
            #print('duplex_rate:', duplex_rate)

    return simplex_reads_basecalled, simplex_reads_filtered, duplex_reads_basecalled, duplex_rate

# Function to create the plots
def create_plots(simplex_reads_basecalled, simplex_reads_filtered, duplex_reads_basecalled, duplex_rate, output_file):
    # Calculate totals
    total_reads = simplex_reads_basecalled + simplex_reads_filtered + duplex_reads_basecalled

    # Create plots
    fig, axs = plt.subplots(1, 2, figsize=(6, 5))

    # Left plot: Duplex rate percentage
    axs[0].bar(1, duplex_rate, color='lightblue', label='Duplex')
    axs[0].bar(1, 100 - duplex_rate, bottom=duplex_rate, color='lightcoral', label='Simplex')
    axs[0].set_ylim(0, 100)
    axs[0].set_ylabel('Percentage')
    axs[0].set_title('Duplex Rate Percentage')
    axs[0].set_xticks([1])
    axs[0].set_xticklabels(['Duplex Rate'])
    axs[0].legend()

    # Right plot: Reads number
    axs[1].bar(1, duplex_reads_basecalled, color='lightblue', label='Duplex')
    axs[1].bar(1, simplex_reads_basecalled + simplex_reads_filtered, bottom=duplex_reads_basecalled, color='lightcoral', label='Simplex')
    axs[1].set_ylim(0, total_reads)
    axs[1].set_ylabel('Number of Reads')
    axs[1].set_title('Reads Number')
    axs[1].set_xticks([1])
    axs[1].set_xticklabels(['Reads'])
    axs[1].legend()

    # Adjust y-axis ticks to show actual numbers
    axs[1].yaxis.set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))

    plt.tight_layout()
    plt.savefig(output_file)
    plt.show()

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Process log file for duplex and simplex read data.')
    parser.add_argument('logfile', type=argparse.FileType('r'), help='Path to the log file')

    # Parse arguments
    args = parser.parse_args()

    # Read log file content
    log_content = args.logfile.read()

    # Extract data from log content
    simplex_reads_basecalled, simplex_reads_filtered, duplex_reads_basecalled, duplex_rate = extract_data(log_content)

    # Create output filename
    log_filename = os.path.basename(args.logfile.name)
    output_file = os.path.splitext(log_filename)[0] + '_duplex_rate.png'

    # Create plots
    create_plots(simplex_reads_basecalled, simplex_reads_filtered, duplex_reads_basecalled, duplex_rate, output_file)

if __name__ == '__main__':
    main()
