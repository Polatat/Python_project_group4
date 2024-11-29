import argparse
import logging
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
from src.parsing.parsing_fastq import FastqParser
from src.statistic.statistic import FastqStat
from src.filter.filter  import filter_fastq

#Set up logging
def setup_logging(output_dir: str) -> None:

    log_file = os.path.join(output_dir,'processing.log')
    logging.basicConfig(
        level= logging.INFO,
        format= '%(asctime)s - %(levelname)s - %(message)s',
        handlers =[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )
    logging.info("Logging has been configured.")

#Set Arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description='Filter and analyze FASTQ data')

    
    parser.add_argument(
                        '-i', '--input', 
                        required =True , 
                        metavar= 'INPUT_PATH',
                        help =' Please choose path to  your input FASTQ file'
                        )
    
    parser.add_argument( 
                        '-o', '--output_dir',
                        required=True,
                        metavar= 'OUTPUT_PATH',
                        help= 'Please choose or name specific directory to save your output files'
                        )

    parser.add_argument(
                        '-qc', '--quality_threshold',
                        metavar= 'QUALITY_SCORE_THRESHOLD',
                        type=float,
                        default=20,
                        help = 'The minimum average quality score to retain a sequence. Default is 20'
                        )      

    parser.add_argument(
                        '-l' , '--min_length',
                        metavar= 'MINIMUM_LENGTH',
                        type= int,
                        default=50,
                        help ='The minimum sequence length to retain. Default is 50'
                       )

    parser.add_argument(
                        '-gc_min', '--gc_minimum',
                        type= float,
                        default = 30.0,
                        help =' The minimum GC content percentage to reatain a sequence. Default is 30.0'
                        )

    parser.add_argument (
                        '-gc_max','--gc_maximum',
                        type = float,
                        default =60.0,
                        help = 'The maximum GC content percentage to retain a sequence. Default is 60.0'
                      )

    parser.add_argument (
                        '-bc_l' , '--barcode_length',
                        metavar= 'LENGTH_OF_BARCODE_SEQUENCE',
                        type=int,
                        default=6,
                        help='The length of the barcode sequence. Default is 6'
                        )

    parser.add_argument(
                        '--header_barcode',
                        action='store_true',
                        default =False,
                        help =' Enable this function to Identify  group of barcode at header sequence'
                        )
    return parser.parse_args()                    

# Pie chart
def pie_chart (distribution, title, output):   #distribution = dis 
    labels = list(distribution.keys())
    sizes = list(distribution.values())

    colors =plt.cm.Paired(range(len(labels)))

    plt.figure(figsize =(8,8))
    plt.pie(sizes, 
            labels=labels, 
            colors =colors, 
            autopct ='%1.2f%%', 
            startangle= 140)

    plt.title(title)
    plt.axis('equal')
    plt.tight_layout()
    plt.savefig(output)
    plt.close()
    logging.info(f"Pie chart is saved to {output}")

#Set up main

def main():

    args = parse_arguments()

    # Checking existed output directory
    os.makedirs(args.output_dir ,exist_ok =True)

    setup_logging(args.output_dir)

    logging.info(f"Beginning processing with FASTQ file : {args.input}")

    # Output file paths
    stat_ori_csv = os.path.join(args.output_dir,'original_statistics.csv')
    stat_fil_csv = os.path.join(args.output_dir,'filtered_statistics.csv')
    piechart_ori = os.path.join(args.output_dir, 'original_barcode_distribtion_piechart.png')
    piechart_fil = os.path.join(args.output_dir,'filtered_barcode_distribution_piechart.png')
    filtered_fastq = os.path.join(args.output_dir,f"filtered.fastq")

    # Step 1: Parsing and Calculating Statistical of Original FASTQ File.
    try:
        logging.info("Beginning FASTQ parsing (Original Data)....")
        parser = FastqParser(args.input)
        parser.parse_fastq()
    
        grouped_sequences = parser.get_grouped_sequences()
        logging.info(f"Parsed {len(grouped_sequences)} groups from the original FASTQ file.")

        logging.info("Generating statistics for Original Data....")
        stat_cal_original = FastqStat(args.input) # args.input = Original FASTq file.
        metrics_original = stat_cal_original.calculate_metrics(
                           barcode_length = args.barcode_length,
                           header_barcode =args.header_barcode
        )
        # Tranform data to CSV file by using pandas
        df_ori = pd.DataFrame(metrics_original)
        df_ori.to_csv(stat_ori_csv, index =False)
        logging.info (f"Original statistics is saved to {stat_ori_csv}.")

        # Create pie chart for orignal barcode distribution
        barcode_dist_original = stat_cal_original.calculate_barcode_distribution(
                                barcode_length= args.barcode_length,
                                header_barcode =args.header_barcode
        )
        pie_chart(
            distribution= barcode_dist_original,
            title       = 'Original Barcode Distribution',
            output      = piechart_ori
        )
        
    except Exception as e:
        logging.error(f"Error during parsing or calculating statistics (Original Data) : {e}")
        sys.exit(1)

#Step 2: Filtering
    try:
        logging.info("Beginning FASTQ filtering....")
        total_sequences, passed_sequences = filter_fastq(
            input_file= args.input,
            output_file= filtered_fastq,
            quality_threshold= args.quality_threshold,
            min_length= args.min_length,
            gc_min = args.gc_minimum,
            gc_max= args.gc_maximum
        )
        failed_sequences = total_sequences - passed_sequences 
        passed_percentage = passed_sequences/total_sequences *100 if total_sequences else 0
        failed_percentage = failed_sequences/total_sequences *100 if total_sequences else 0
        logging.info(f'Number of sequences in original FASTQ: {total_sequences}(100.00%)')
        logging.info(f'Number of sequences after filtering: {passed_sequences}({passed_percentage:.2f}%)')
        logging.info(f'Number of sequences failed filtering: {failed_sequences} ({failed_percentage:.2f}%)')

    except Exception as e: 
        logging.error(f" Error during filtering: {e}")
        sys.exit(1)

    # Step 3 :Parsing and Calculating Statistical of Filtered FASTQ File.
    try:
        logging.info("Beginning FASTQ parsing (Filtereed Data)....")
        parser_filtered = FastqParser(filtered_fastq)
        parser_filtered.parse_fastq()
        grouped_sequences_filtered = parser_filtered.get_grouped_sequences()
        logging.info( f"Parsed {len(grouped_sequences_filtered)} groups from the filtered FASTQ file.")

        logging.info("Generating statistics for Filtered Data....")
        stat_cal_filtered = FastqStat(filtered_fastq)
        metrics_filtered = stat_cal_filtered.calculate_metrics(
                            barcode_length= args.barcode_length,
                            header_barcode = args.header_barcode
        )

        # Tranform data to CSV file by using pandas
        df_filtered = pd.DataFrame(metrics_filtered)
        df_filtered.to_csv(stat_fil_csv, index= False)
        logging.info(f"Filtered statistics saved to {stat_fil_csv}.")

        # Filtered barcode distribution pie chart
        barcode_dist_filtered = stat_cal_filtered.calculate_barcode_distribution(
                                barcode_length=args.barcode_length,
                                header_barcode=args.header_barcode
        )
        pie_chart(
                distribution= barcode_dist_filtered,
                title ='Filtered Barcode Distribution',
                output = piechart_fil
        )
    except Exception as e:
        logging.error ( f" Error during parsing or calculating statistics (Filtered Data) : {e}")
        sys.exit(1)

    logging.info("Processing Completed Successfully.")

if __name__ == "__main__":
    main()


       
                                        
        





    

















    














    
    

    
    

    

