import os

if __name__ == '__main__':
    base_dir = r"C:\Users\ab50\Downloads"
    file_name = "Cdb.csv"
    bin_dict = {}
    with open(os.path.join(base_dir,file_name),'r') as f_handle:
        f_handle.readline()
        prev_bin_id = '-1'
        for line in f_handle:
            bin_name = line.split(',')[0]
            curr_bin_id = (line.split(',')[1]).split('_')[0]
            if curr_bin_id != prev_bin_id:
                bin_dict[curr_bin_id] = [bin_name]
                prev_bin_id = curr_bin_id
            else:
                bin_dict[curr_bin_id].append(bin_name)
    out_file_name = os.path.join(base_dir,"clustered_bins.csv")
    with open(out_file_name,'w') as f_handle:
        f_handle.write("De-replicated Bin,Cluster Size,Clustered Bin Names\n")
        for key in bin_dict.keys():
            clst_bins = "|".join(bin_dict[key][1:])
            f_handle.write(bin_dict[key][0] + "," + str(len(bin_dict[key])) + "," + clst_bins+"\n")









