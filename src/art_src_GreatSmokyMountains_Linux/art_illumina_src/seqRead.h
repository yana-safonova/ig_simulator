/*
 * >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 * ART -- Artificial Read Transcription, Illumina Q version 
 * Authors: Weichun Huang 2008-2016
 * License: GPL v3 
 * ############################################################################
 * #    This program is free software: you can redistribute it and/or modify  #
 * #    it under the terms of the GNU General Public License as published by  #
 * #    the Free Software Foundation, either version 3 of the License, or     #
 * #    (at your option) any later version.                                   #
 * #                                                                          #
 * #    This program is distributed in the hope that it will be useful,       #
 * #    but WITHOUT ANY WARRANTY; without even the implied warranty of        #
 * #    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         #
 * #    GNU General Public License for more details.                          #
 * #                                                                          #
 * #    You should have received a copy of the GNU General Public License     #
 * #    along with this program.  If not, see <http://www.gnu.org/licenses/>. #
 * ############################################################################
 * <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*/

#pragma once

#include <cmath>
#include <ctime>
#include <string>
#include <iostream>
#include <vector>
#include <iterator>
#include <map>
#include <cstdio>
#include <set>
#include <gsl/gsl_cdf.h>

#include "empdist.h"

using namespace std;

class seqRead{
public:
    seqRead(){  }

    static const int max_N_qual=10;
    vector <double> ins_rate; //Bionomial cum_prob gsl_cdf_binomial_Q(unsigned int k, double p, unsigned int n) 
    vector <double> del_rate; //Binomial
    vector<double> sub_rate; //Binomial
    void set_rate(int read_len, double p, int max_num, vector <double>& rate){
        rate.resize(max_num);
        for(size_t i=1; i<=max_num; i++){
            rate[i-1]= gsl_cdf_binomial_Q(i, p, read_len);
        }
    };
    static char rand_base(){
        short base=(short)ceil(r_prob()*4);
        switch(base){
            case 1:
                return 'A';
            case 2:
                return 'C';
            case 3:
                return 'G';
            case 4:
                return 'T';  
	    default:
                return 'N';  
        }
    };

    //static bool with_indel;
    int get_indel(int read_len);
    //number of deletions <= number of insertions
    int get_indel_2(int read_len);
    map<int,char,less<int> > indel;
    map<int,char> substitution;
    bool is_plus_strand;
    unsigned long bpos; //parent
    string seq_read;
    string seq_ref;
    void clear(){
        indel.clear();
        substitution.clear();
        seq_read.clear();
        seq_ref.clear();
    }
    //string aln_read;
    //string aln_ref;
    bool get_aln(string& aln_read, string& aln_ref){
        if(indel.size()==0) return false;
        map<int,char,less<int> >::iterator it;
        aln_read=seq_read; aln_ref=seq_ref;
        for(it=indel.begin(); it!=indel.end(); it++){
            if(it->second!='-'){
                aln_ref.insert(it->first,1,'-');
            }
            else{
                aln_read.insert(it->first,1,'-');
            }
        }
        return true;
    };

    void ref2read();

    //based on empirical dist of quali scores
    int add_error(vector<short>&qual){
        if(qual.size()!=seq_read.size()){
            cerr<<"Error: the number of bases is not equal to the number of quality scores!\n";
	    cerr<< "qual size: " << qual.size() << ",  read len: " << seq_read.size() << endl;
	    exit(1);
            return 0;
        }
        int num=0;
        for(size_t i=0; i<qual.size(); i++){
	    if(seq_read[i]=='N'){
		    //qual[i]=(short) floor(max_N_qual*r_prob()); //reassign N base quality score < max_N_qual   
		    qual[i]=(short)1;
		    continue;
	    }
            if(r_prob()<empdist::prob_err[qual[i]]){
                char achar=seq_read[i];
                while(seq_read[i]==achar){ achar=rand_base(); }
                seq_read[i]=achar;
                substitution[i]=achar;
                num++;
            }
        }
        return num;
    };

    //based on empirical dist of quali scores
    int add_error_stat(vector<short>&qual, vector<long>&err_pos){
        if(qual.size()!=seq_read.size()){
            cerr<<qual.size()<<"\t"<<seq_read.size()<<endl;
            cerr<<"Error: the number of bases is not equal to the number of quality score!\n";
	    exit(1);
        }
        int num=0;
        for(size_t i=0; i<qual.size(); i++){
            if(r_prob()<empdist::prob_err[qual[i]]){
              err_pos[i]+=1;
                char achar=seq_read[i];
                while(seq_read[i]==achar){ achar=rand_base(); }
                seq_read[i]=achar;
                substitution[i]=achar;
                num++;
            }
        }
        return num;
    };

    //base on bionomial error substitution rate
    int add_error(int read_len){
        int sub_num=0;
        for(int i=(int)sub_rate.size()-1; i>=0; i--){
            if(sub_rate[i]>=r_prob()){
                sub_num=i+1;
                for(int j=i; j>=0;){
                    int pos=(int) floor(read_len*r_prob());
                    if(substitution.count(pos)==0){
                        char achar=seq_read[pos];
                        while(seq_read[pos]==achar){ achar=rand_base(); }
                        substitution[pos]=achar;
                        seq_read[pos]=achar;
                        j--;
                    }
                }
                break;
            }
        }
        return sub_num;
    };

    //vector<char> aln_read;
    //vector<char> aln_ref;
};
