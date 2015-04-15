#include "../utils/include_me.hpp"
#include "../utils/fastq_reader.hpp"
#include "../utils/cluster_utils/clusterizator.hpp"
#include "../utils/string_tools.hpp"
#include "omp.h"
#include "../utils/cluster_utils/clusters.hpp"

string GetResultFname(string reads_fname) {
	size_t pos = reads_fname.find(".fastq");
	string base_name = reads_fname.substr(0, reads_fname.size() - 6);
	if(pos == string::npos) {
		pos = reads_fname.find(".fq");
		base_name = reads_fname.substr(0, reads_fname.size() - 3);
	}
	return base_name;
}

int main(int argc, char *argv[]) {

    if(argc != 3) {
        cout << "Ivalid input paramters" <<
            "\targv[1] - FASTQ with merged simulated reads" << endl <<
            "\targv[2] - basename for output files" << endl;
        return 1;
    }

	string fname(argv[1]);
	vector<fastq_read> reads = SingleFastqReader(fname).ReadFile();
	cout << reads.size() << " reads were extracted from " << argv[1] << endl;

	Clusterization clusters(reads);
	string prefix = "ID=";
    map<size_t, pair<string, size_t> > cluster_seq_size;
	for(size_t i = 0; i < reads.size(); i++) {
			vector<string> splits = split(reads[i].name, '_');
			string clust_num = splits[4];
//			splits = split(clust_num, '_');
//			size_t prefix_pos = splits[3].find(prefix);
//			assert(prefix_pos != string::npos);
//			clust_num = splits[3].substr(prefix_pos + prefix.size(), clust_num.size() - prefix.size() - prefix_pos);
            size_t cluster_id = string_to_number<size_t>(clust_num);
			clusters.Add(i, cluster_id);
            
            if(cluster_seq_size.find(cluster_id) == cluster_seq_size.end()) {
                cluster_seq_size[cluster_id].first = reads[i].seq;
                cluster_seq_size[cluster_id].second = 1;
            }
            else
                cluster_seq_size[cluster_id].second++;
		}
    string basename(argv[2]);
	clusters.WriteReadClusterMap(basename + ".rcm");
	//cout << "Read-cluster map was written to " << basename + ".rcm" << endl;

    // filling cluster FASTA information
    for(auto it = cluster_seq_size.begin(); it != cluster_seq_size.end(); it++)
        clusters.AddClusterSequence(it->first, it->second.first, it->second.second);

	clusters.WriteClustersFasta(basename + ".clusters.fa");
	//cout << "Read-cluster map was written to " << basename + ".clusters.fa" << endl;
}
