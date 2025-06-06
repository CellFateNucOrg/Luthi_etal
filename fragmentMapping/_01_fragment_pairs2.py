from collections import defaultdict
#from threading import Thread, Lock
#from itertools import count
import argparse
import bioframe as bf


def count_hic_pairs(filename, output_filename, min_distance, max_distance, regions=None):
  """
  Counts occurrences of HIC pairs in a large file using multiple threads.

  Args:
    filename: Path to the file containing HIC data.
    output_filename: Path to output file
  """
  print("parsing contacts from file: "+filename)
  hic_pair_counts = defaultdict(int)
  with open(filename, 'r') as file:
    try:
        line = file.readline()
        while line:
            # Split line by tabs
            parts = line.strip().split('\t')
            chr1=parts[1]
            start1=int(parts[2])
            chr2=parts[4]
            start2=int(parts[5])
            end2=int(parts[5])+1
            df2 = bf.from_any([[chr2, start2, end2]], name_col='chrom')
            # Extract and sort HIC regions (avoiding unnecessary splits)
            #print("Pair 1: " + parts[8])
            hic1 = parts[8]
            #print(hic1)
            #print("Pair 2: " + parts[9])
            hic2 = parts[9]
            hic_pair = min(hic1, hic2) + '|' + max(hic1, hic2)
            #print(hic_pair)
            if chr1==chr2 and abs(start2-start1) > min_distance and abs(start2-start1) < max_distance: # only if in cis
                if regions is not None:
                  end1=int(parts[2])+1
                  df1 = bf.from_any([[chr1, start1, end1]], name_col='chrom')
                  end2=int(parts[5])+1
                  df2 = bf.from_any([[chr2, start2, end2]], name_col='chrom')
                  anchor1overlap = len(bf.overlap(df1,regions, how="inner"))  
                  anchor2overlap = len(bf.overlap(df2, regions, how="inner"))
                  if anchor1overlap > 0 and anchor2overlap > 0:
                    hic_pair_counts[hic_pair] += 1
                else:
                  hic_pair_counts[hic_pair] += 1
            line=file.readline()
    except Exception as e:
        print(f"Error processing line: {e}")


  # Save HIC pairs and their counts to a file (can be modified as needed)
  # ... (same logic as before)
    # Save HIC pairs and their counts to a file
  with open(output_filename, 'w') as output_file:
    for hic_pair, count in hic_pair_counts.items():
      output_file.write(f"{hic_pair} {count}\n")
  

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='Extract valid pairs')
  parser.add_argument('-i', '--input_filename')
  parser.add_argument('-o', '--output_filename')
  parser.add_argument('-r', '--regions', default=None)
  parser.add_argument('-n', '--min_distance', type=int, default=0) 
  parser.add_argument('-d', '--max_distance', type=int, default=30_000)
  args = parser.parse_args()
  input_filename = args.input_filename
  output_filename = args.output_filename
  min_distance = args.min_distance
  max_distance = args.max_distance
  if args.regions is not None:
    regions = bf.read_table(args.regions,schema="bed6")
  else:
    regions = None
  count_hic_pairs(input_filename, output_filename, min_distance, max_distance, regions)

