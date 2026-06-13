import os
import shutil

out = os.system('export PGAP_INPUT_DIR=/scratch/gpfs/DONIA/abiswas/tools/pgap/.pgap')

if __name__ == '__main__':
    pfam_name = 'AAC6__III-homologs'
    data_dir = '/scratch/gpfs/DONIA/abiswas/MetaBGCRuns/' + pfam_name + '/output/analytics_10_10/scf_blast_results/scf_match' 
    script_dir = '/scratch/gpfs/DONIA/abiswas/annotations/' + pfam_name + '/scripts'
    pgap_input_dir = '/scratch/gpfs/DONIA/abiswas/annotations/' + pfam_name + '/run_folders' 
    pgap_output_dir = '/scratch/gpfs/DONIA/abiswas/annotations/' + pfam_name + '/annotations'

    os.makedirs(pfam_name, exist_ok=True)
    os.makedirs(pgap_input_dir, exist_ok=True)
    os.makedirs(pgap_output_dir, exist_ok=True)
    os.makedirs(script_dir, exist_ok=True)
    os.makedirs(os.path.join(script_dir, 'logs'), exist_ok=True)
    shutil.copy2('submol.yaml', script_dir)
    shutil.copy2('runPGAPArray.sh', script_dir)
    shutil.copy2('runPGAP.sh', script_dir)

    with open(os.path.join(script_dir, 'pgap_all_samples') , 'w') as  sample_run_file:
        for sample_file in os.listdir(data_dir):
            f_path = os.path.join(data_dir, sample_file)
            sample_name = os.path.splitext(sample_file)[0]
            # Check scaffold file is not empty
            if not os.path.getsize(f_path) == 0:
                # Check if output directory exists and run is completed
                # Remove the old old directory if there but no annotation
                pgap_sample_out_dir = os.path.join(pgap_output_dir, sample_name + '_results')
                if os.path.exists(pgap_sample_out_dir) and os.path.isdir(pgap_sample_out_dir):
                    faa_file = os.path.join(pgap_sample_out_dir, 'annot.faa')
                    if not (os.path.exists(faa_file) and os.path.getsize(faa_file) > 0):
                        shutil.rmtree(pgap_sample_out_dir)
                    else:
                        continue
                # Add to file to run
                sample_run_file.write(sample_name + '\n')
                # Create input files 
                run_dir = os.path.join(pgap_input_dir, sample_name)
                os.makedirs(run_dir, exist_ok=True)
                shutil.copy2(f_path, os.path.join(run_dir, sample_file))
                shutil.copy2(os.path.join(script_dir, 'submol.yaml'), os.path.join(run_dir, 'submol.yaml'))
                with open(os.path.join(run_dir, 'input.yaml'),'w') as f_yml:
                    f_yml.write('fasta:\n  class: File\n  location: ' + sample_file  + '\n')
                    f_yml.write('submol:\n  class: File\n  location: submol.yaml\n')

