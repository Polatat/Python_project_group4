   
import sys
from Bio import SeqIO
import logging

def calculate_gc_content(seq):
    gc = seq.count('G') + seq.count('C') + seq.count('g') + seq.count('c')
    return (gc / len(seq)) * 100 if len(seq) > 0 else 0

def filter_fastq(input_file, output_file,quality_threshold=20, min_length=50, gc_min=30, gc_max=60):
    total = 0
    passed = 0
    
    if quality_threshold is None: 
        logging.info( f" Quality threshold  not detected. Default value is 20")
        quality_threshold = 20
    else:
        logging.info(f"Quality threshold set to {quality_threshold}")

    with open(input_file, "r") as in_handle, open(output_file, "w") as out_handle:
        for record in SeqIO.parse(in_handle, "fastq"):
            total += 1
            # Calculate average quality score
            qualities = record.letter_annotations.get("phred_quality",[])
            if not qualities:
                logging.debug(f" Sequence{record.id} has no quality scores. Skipping.")
                continue
            avg_quality = sum(qualities)/ len(qualities)

            if avg_quality < quality_threshold:
                continue
            
            # Check sequence length
            seq_length = len(record.seq)
            if seq_length <= min_length:
                continue

            # Check GC content
            gc_content = calculate_gc_content(str(record.seq))
            if gc_content < gc_min or gc_content > gc_max:
                continue
            
            # Write to output if all conditions are met
            SeqIO.write(record, out_handle, "fastq")
            passed += 1
            
            # Log progress every 100,000 sequences
            # if total % 100000 == 0:
            #     logging.info(f"Processed {total} sequences, {passed} passed so far.")
    
    logging.info(f"Filtering complete. Total sequences: {total}, Passed: {passed}, Failed: {total - passed}")
    return total, passed

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python filter_fastq.py input.fastq output_filtered.fastq")
        sys.exit(1)
    
    input_fastq = sys.argv[1]
    output_fastq = sys.argv[2]
    
    # Set up logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s')
    
  # Customize criteria
    QUALITY_THRESHOLD = 20  # or 30
    MIN_LENGTH = 50
    GC_MIN = 30            # Minimum GC content percentage
    GC_MAX = 60             # Maximum GC content percentage
    
    logging.info("Starting FASTQ filtering...")
 
    total_sequences, passed_sequences = filter_fastq(
        input_fastq, 
        output_fastq, 
       quality_threshold= QUALITY_THRESHOLD, 
       min_length= MIN_LENGTH, 
       gc_min=GC_MIN, 
       gc_max=GC_MAX
    )
    

    # failed_sequences = total_sequences - passed_sequences 
    # total_percentange =  total_sequences/total_sequences * 100
    # passed_percentage = passed_sequences/total_sequences *100
    # failed_percentage = failed_sequences/total_sequences *100 
    # logging.info(f'Number of sequences in original FASTQ: {total_sequences}({total_percentange:.2f}%)')
    # logging.info(f'Number of sequences after QC: {passed_sequences}({passed_percentage:.2f}%)')
    # logging.info(f'Number of sequences failed QC: {failed_sequences} ({failed_percentage:.2f}%)')
    




    
