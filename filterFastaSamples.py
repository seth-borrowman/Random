import argparse

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('-i', '--input', help='input FASTA', required=True)
  parser.add_argument('-l', '--list', help='File with a record name to keep on each line', required=True)
  parser.add_argument('-o', '--output', help='output FASTA', required=False, default='filtered.fasta')
  args = parser.parse_args()

  in_fasta=open(args.input)
  out_fasta=open(args.output, 'w')
  sample_list=open(args.list)

  samples_to_keep = sample_list.readlines()
  print(samples_to_keep)
  sample_list.close()
  in_fasta_lines = in_fasta.readlines()

  keep_record = False
  for i in in_fasta_lines:
    if i.startswith('>'):
      if i[1:] in samples_to_keep:
        keep_record = True
        out_fasta.write(i)
      else:
        keep_record = False
    else:
      if bool(keep_record):
        out_fasta.write(i)
  
  out_fasta.close()
  in_fasta.close()

if __name__ == '__main__':
  main()
