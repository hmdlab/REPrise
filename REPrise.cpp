/*
REPrise_v1.0.1.cpp
Author: Atsushi Takeda, Daisuke Nonaka
*/


#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <queue>
#include <set>
#include <string>
#include <sys/resource.h>
#include <sys/time.h>
#include <tuple>
#include <unordered_map>
#include <vector>

using namespace std;

#include "cmd_line_opts.h"
#include "sais_long.h"
#include "REPrise.hpp"

#define IUPAC(c) c == 'R' || c == 'r' || c == 'Y' || c == 'y' || c == 'M' || c == 'm' || c == 'K' || c == 'k' || c == 'W' || c == 'w' || c == 'S' || c == 's' || c == 'B' || c == 'b' || c == 'D' || c == 'd' || c == 'H' || c == 'h' || c == 'V' || c == 'v'

string SEQUENCE_FILE;
string FREQ_FILE;
string REPROF_FILE;
string MASKED_FILE;
string BED_FILE;//20211201追加

#define seq_type long //20210701修正,pos修正

//

int TANDEMDIST;
int MINFREQ;
int MINIMPROVEMENT;
double MAXENTROPY = -0.7;
int PADLENGTH = 11000; /* should be >= L+OFFSETWIDTH */
int MAXEXTEND;         /* How far to extend. max total length of consensus is 2*MAXEXTEND+l  (MAXEXTEND >>> l) */
int MAXREPEAT;
int OFFSETWIDTH;  /* max offset (5) */
int WHEN_TO_STOP; /* stop if no improvement after extending this far (500) */
int MINLENGTH;
int CAPPENALTY;
int MATCHSCORE;
int MISMATCHSCORE;
int GAPSCORE;
int GAPEXTENDSCORE;

int MINSCORE = -100000;

int VERBOSE;
int OUTPUT_MASKED_REGION;
int HELP;

char* consensus;
//seq_type* pos;
//bool* rev;
vector<seq_type> pos;
vector<bool> rev;

seq_type length;
int k;
char* sequence;
seq_type* SA;

vector<pair<string, seq_type> > chrtable;

int PARALLELNUM;

#ifdef _OPENMP
double time1, totaltime;
#else
clock_t time1, totaltime;
#endif

int main(int argc, char* argv[]) {
#ifdef _OPENMP
    time1 = omp_get_wtime();
#else
    time1 = clock();
#endif
    totaltime = 0;

    int KMERDIST;
    char* input;
    char* output;

    co_get_int(argc, argv, "-match", &MATCHSCORE) || (MATCHSCORE = 1);
    co_get_int(argc, argv, "-mismatch", &MISMATCHSCORE) || (MISMATCHSCORE = -1);
    co_get_int(argc, argv, "-gap", &GAPSCORE) || (GAPSCORE = -5);
    co_get_int(argc, argv, "-gapex", &GAPEXTENDSCORE) || (GAPEXTENDSCORE = -1);
    co_get_int(argc, argv, "-cappenalty", &CAPPENALTY) || (CAPPENALTY = -20);
    co_get_int(argc, argv, "-dist", &KMERDIST) || (KMERDIST = 0);

    co_get_int(argc, argv, "-maxextend", &MAXEXTEND) || (MAXEXTEND = 10000);
    co_get_int(argc, argv, "-maxrepeat", &MAXREPEAT) || (MAXREPEAT = 100000);
    co_get_int(argc, argv, "-maxgap", &OFFSETWIDTH) || (OFFSETWIDTH = 5);
    co_get_int(argc, argv, "-stopafter", &WHEN_TO_STOP) || (WHEN_TO_STOP = 100);
    co_get_int(argc, argv, "-minlength", &MINLENGTH) || (MINLENGTH = 50);
    co_get_int(argc, argv, "-minfreq", &MINFREQ) || (MINFREQ = 3);
    co_get_int(argc, argv, "-minimprovement", &MINIMPROVEMENT) || (MINIMPROVEMENT = 3);
    co_get_int(argc, argv, "-tandemdist", &TANDEMDIST) || (TANDEMDIST = 500);

    co_get_int(argc, argv, "-pa", &PARALLELNUM) || (PARALLELNUM = 1);

    co_get_bool(argc, argv, "-verbose", &VERBOSE) || (VERBOSE = 0);
    co_get_bool(argc, argv, "-h", &HELP) || (HELP = 0);
    co_get_bool(argc, argv, "-additonalfile", &OUTPUT_MASKED_REGION) || (OUTPUT_MASKED_REGION = 0);

    if(HELP==1){
        print_usage();
        exit(1);
    }

    co_get_string(argc, argv, "-input", &input);
    co_get_string(argc, argv, "-output", &output);
    SEQUENCE_FILE = input;
    FREQ_FILE = output;
    FREQ_FILE += ".freq";
    REPROF_FILE = output;
    REPROF_FILE += ".reprof";
    MASKED_FILE = output;
    MASKED_FILE += ".masked";
    BED_FILE = output;
    BED_FILE += ".bed";

    build_sequence();
    //printf("%d\n", VERBOSE);
    //debug for chrtrace
    for(int i = 0; i < chrtable.size(); i++){
        cout << chrtable[i].first << "\t" << chrtable[i].second << endl;
    }

    co_get_int(argc, argv, "-k", &k) || (k = default_k(length,KMERDIST));
    cout << "kmer length: " << k << endl;

    SA = (seq_type*)malloc(length * sizeof(seq_type));
    sais(sequence, SA, length);

    int dist0cachelength = 1 << min(min(k, default_k(length, KMERDIST) - 2), 13) * 2;
    int cachelength = 1 << min(min(k, default_k(length,KMERDIST) - 2), 13 - KMERDIST * 2) * 2;
    /*20210415変更*/
    vector<vector<tuple<seq_type, seq_type, char, char>>> dist0cachetable(dist0cachelength);
    vector<vector<tuple<seq_type, seq_type, char, char>>> cachetable(cachelength); 

    store_cache(0, dist0cachetable);
    display_time("dist0cache");
    if (KMERDIST > 0) {
        store_cache(KMERDIST, cachetable);
        display_time("nomalcache");
    }

    priority_queue<pair<int, vector<char>>> kmers;
    if (KMERDIST > 0)
        build_sortedkmers(kmers, dist0cachetable, cachetable);
    else
        build_sortedkmers(kmers, dist0cachetable, dist0cachetable);
    if (KMERDIST > 0)
        build_repeat_families(kmers, cachetable);
    else
        build_repeat_families(kmers, dist0cachetable);

    return 0;
}

void store_cache(int edit_distance, vector<vector<tuple<seq_type, seq_type, char, char>>>& cachetable) {
    seq_type cachelength = 0;
    seq_type sz = cachetable.size();
    while (sz > 3) {
        sz = sz >> 2;
        cachelength++;
    }
    sz = cachetable.size();

#pragma omp parallel for
    for (int i = 0; i < PARALLELNUM; i++) {
        vector<char> query(cachelength);
        set<tuple<seq_type, seq_type, char, char>> matched;
        for (long j = sz * i / PARALLELNUM; j < sz * (i + 1) / PARALLELNUM; j++) {
            for (int l = cachelength - 1; l >= 0; l--)
                query[l] = (j >> l * 2) & 3;
            SA_search(query, 0, length, 0, 0, edit_distance, matched);
            cachetable[j].clear();
            for (const auto& e : matched)
                cachetable[j].push_back(e);
            matched.clear();
        }
    }
}

class ComparePair {
public:
    bool operator()(pair<seq_type, seq_type> n1, pair<seq_type, seq_type> n2) {
        if (n1.first != n2.first)
            return n1.first < n2.first;
        else
            return n1.second >= n2.second;
    }
};

void build_sortedkmers(priority_queue<pair<int, vector<char>>>& kmers, const vector<vector<tuple<seq_type, seq_type, char, char>>>& dist0cachetable, const vector<vector<tuple<seq_type, seq_type, char, char>>>& cachetable) {
    ofstream fout;
    fout.open(FREQ_FILE, ios::out);
    if (!fout) {
        cerr << "Could not open " << FREQ_FILE << endl;
        exit(1);
    }

    vector<bool> samekmer_flag(length, false);
    vector<bool> invalidkmer_flag(length, false);
    vector<bool> entropy_flag(length, false);
    vector<bool> mask_flag(length, false);

    vector<char> query(k);
    vector<seq_type> occs;
    for (seq_type i = 0; i < length; i++) {
        if (i > length - k)
            invalidkmer_flag[i] = true;
        else {
            for (int j = 0; j < k; j++) {
                query[j] = sequence[j + i];
                if (query[j] > 3)
                    invalidkmer_flag[i] = true;
            }
        }
        if (invalidkmer_flag[i])
            continue;
        if (compute_entropy(query) > MAXENTROPY) {
            entropy_flag[i] = true;
            continue;
        }
        if (!samekmer_flag[i]) {
            occs = findkmer(query, dist0cachetable);
            for (const auto& e : occs)
                samekmer_flag[e] = true;
            occs = findkmer(reverse_complement(query), dist0cachetable);
            for (const auto& e : occs)
                samekmer_flag[e] = true;
            samekmer_flag[i] = false;
        }
    }
    display_time("buildflags");

    vector<priority_queue<pair<int, seq_type>>> freqocc(PARALLELNUM);

    //int loopnum = max(1, length / 100000);
#pragma omp parallel for
    for (int i = 0; i < PARALLELNUM; i++) {
        //for (int h = 0; h < loopnum; h++) {
        vector<char> query(k);
        vector<seq_type> occs;
        int freq;
        for (seq_type j = ((long long int)length /* h + (long long int)length*/* i / PARALLELNUM) /*/ loopnum*/; j < ((long long int)length /* h + (long long int)length */* (i + 1) / PARALLELNUM) /*/ loopnum*/; j++) {
            if (samekmer_flag[j] || invalidkmer_flag[j] || entropy_flag[j])
                continue;
            for (int l = 0; l < k; l++)
                query[l] = sequence[j + l];
            occs = findkmer(query, cachetable);
            removetandem(occs);
            freq = occs.size();
            occs = findkmer(reverse_complement(query), cachetable);
            removetandem(occs);
            freq += occs.size();
            if (freq >= MINFREQ)
                freqocc[i].push(make_pair(freq, j));
            //}
        }
    }

    display_time("searchkmer");

    priority_queue<pair<seq_type, seq_type>, vector<pair<seq_type, seq_type>>, ComparePair> allfreqocc;
    for (int i = 0; i < PARALLELNUM; i++) {
        while (!freqocc[i].empty()) {
            allfreqocc.push(freqocc[i].top());
            freqocc[i].pop();
        }
    }

    display_time("kmersmerge");
    if (allfreqocc.size() == 0) {
        cout << "no kmer detected" << endl;
        exit(0);
        //if we do not exit here bug in next phase
    }
    vector<seq_type> sortedoccs(allfreqocc.size());
    for (auto& e : sortedoccs) {
        e = allfreqocc.top().second;
        allfreqocc.pop();
    }
    /*以下謎!!!*/
    int BLOCK = min((int)sortedoccs.size(), PARALLELNUM * 1000);//BLOCK化による高速化
    int freq;
    vector<seq_type> tmpoccs(BLOCK);
    vector<vector<seq_type>> occsvec(BLOCK);
    vector<vector<seq_type>> rcoccsvec(BLOCK);
    for (int i = 0; i < (int)sortedoccs.size() / BLOCK; i++) {
#pragma omp parallel for
        for (int j = 0; j < PARALLELNUM; j++) {
            for (int l = 0; l < BLOCK / PARALLELNUM; l++) {
                int index = j + l * PARALLELNUM;
                tmpoccs[index] = sortedoccs[i * BLOCK + index];
                vector<char> query(k);
                for (int m = 0; m < k; m++)//2回目の計算??
                    query[m] = sequence[tmpoccs[index] + m];
                occsvec[index] = findkmer(query, cachetable);
                removetandem(occsvec[index]);
                rcoccsvec[index] = findkmer(reverse_complement(query), cachetable);
                removetandem(rcoccsvec[index]);
            }
        }
        for (int j = 0; j < BLOCK; j++) {
            removemasked(occsvec[j], mask_flag, false);
            removemasked(rcoccsvec[j], mask_flag, true);
            maskbyseed(occsvec[j], mask_flag, false);
            maskbyseed(rcoccsvec[j], mask_flag, true);
            freq = occsvec[j].size() + rcoccsvec[j].size();
            if (freq >= MINFREQ)
                allfreqocc.push(make_pair(freq, tmpoccs[j]));
        }
    }

    display_time("maskbykmer");
    seq_type occ;
    while (!allfreqocc.empty()) {//kmerのファイルへの書き込み
        occ = allfreqocc.top().second;
        freq = allfreqocc.top().first;
        for (int i = 0; i < k; i++)
            query[i] = sequence[occ + i];
        kmers.push(make_pair(freq, query));
        allfreqocc.pop();
        for (int i = 0; i < k; i++)
            fout << num_to_char(sequence[occ + i]);
        fout << "\t" << freq << "\t" << occ << endl;
    }
    fout.close();
    display_time("outputkmer");
}

void build_repeat_families(priority_queue<pair<int, vector<char>>>& kmers, const vector<vector<tuple<seq_type, seq_type, char, char>>>& cachetable) {
    ofstream fout1, fout2, fout3;
    fout1.open(REPROF_FILE, ios::out);
    if (!fout1) {
        cerr << "Could not open " << REPROF_FILE << endl;
        exit(1);
    }
    if(OUTPUT_MASKED_REGION == 1){
        fout2.open(MASKED_FILE, ios::out);
        if (!fout2) {
            cerr << "Could not open " << MASKED_FILE << endl;
            exit(1);
        }
        fout3.open(BED_FILE, ios::out);
        if (!fout3) {
            cerr << "Could not open " << BED_FILE << endl;
            exit(1);
        }
    }


    vector<bool> mask_flag(length, false);
    vector<bool> mask_flag2(length, false);
    allocate_space();
    pair<int, vector<char>> currentseed;
    int seedfreq;
    int consensusstart, consensusend;
    pair<string, seq_type> chr_start;//20211205_takeda_edit gff出力用, 並列化の時はループの内側で宣言する必要があるかも
    int repeat_num = 0;
    int element_count;
    while (repeat_num < MAXREPEAT) {
        element_count = 0;
        if (kmers.empty())
            break;
        else
            currentseed = find_bestseed(kmers, cachetable, mask_flag);
        seedfreq = currentseed.first;
        if (seedfreq < MINFREQ)
            break;

        element_count = 0;
        vector<seq_type> repeatstart(seedfreq);
        vector<seq_type> repeatend(seedfreq);

#pragma omp parallel sections num_threads(2)
        {
#pragma omp section
            {
                vector<seq_type> seed_ext(seedfreq, -1);
                consensusend = MAXEXTEND + extend(true, seedfreq, seed_ext);
                for (int i = 0; i < seedfreq; i++)
                    repeatend[i] = seed_ext[i];
            }
#pragma omp section
            {
                vector<seq_type> seed_ext(seedfreq, -1);
                consensusstart = MAXEXTEND - 1 - extend(false, seedfreq, seed_ext);
                for (int i = 0; i < seedfreq; i++)
                    repeatstart[i] = -seed_ext[i] - 1;
            }
        }


        //20211105_takeda_edit, paiwwise REalignment between consensus-repeat and canditates 
        pair<seq_type, seq_type> repeat_element;//アラインメントの結果を管理する, 並列化の時はループの内側で宣言する必要があるかも...?
        if(consensusend - consensusstart + 1 >= MINLENGTH){
            /*mask_flag, mask_flag2の両方を更新, consensusの出力もする*/
            /*並列化について検討*/
            
            if(VERBOSE == true){
                cout << "-------------------------------------------------------------------------------"<< endl;
                cout << ">R=" << repeat_num << " seedfreq=" << seedfreq << ", length=" << consensusend - consensusstart + 1 << ", Seed=";
                for (int i = 0; i < k; i++)
                    cout << num_to_char(currentseed.second[i]);
                cout << endl;
                seq_type x;
                for (x = consensusstart; x <= consensusend; x++) {
                    cout << num_to_char(consensus[x]);
                    if ((x - consensusstart) % 80 == 79)
                        cout << endl;
                    }
                if ((x - consensusstart) % 80 > 0)
                    cout << endl;
                cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"<< endl;
            }
            
            
            for(int i = 0; i < seedfreq; i++){
                repeat_element = masking_align(i, consensusstart, consensusend);
                //cout << repeat_element.first << "," << repeat_element.second << endl;
                maskbyrepeat_element(i, repeat_element.first, repeat_element.second, mask_flag);
                if(llabs(repeat_element.second - repeat_element.first + 1) < MINLENGTH){//20220212  ,+ 1するように修正
                    continue;
                }else{
                    cout << "!element" << element_count << endl;
                    maskbyrepeat_element(i, repeat_element.first, repeat_element.second, mask_flag2);
                    chr_start = chrtracer(pos[i]);
                    if(OUTPUT_MASKED_REGION == 1){
                        if(!rev[i]){
                            fout3 << chr_start.first << "\t" << pos[i] + repeat_element.first - chr_start.second  << "\t" << pos[i] + repeat_element.second - chr_start.second  << "\t" <<  "R=" << repeat_num  << "\t" << repeat_element.second - repeat_element.first + 1 <<  "\t+" <<  endl; 
                        }else{
                            fout3 << chr_start.first << "\t" << pos[i] + repeat_element.first - chr_start.second  << "\t" << pos[i] + repeat_element.second - chr_start.second  << "\t" <<  "R=" << repeat_num  << "\t" << repeat_element.second - repeat_element.first + 1 <<  "\t-" <<  endl; 
                        }
                    }
                    
                    if(VERBOSE = true){
                            seq_type x;
                    for (x = repeat_element.first; x <= repeat_element.second; x++) {
                        cout << num_to_char(sequence[pos[i] + x]);
                        if ((x - repeat_element.first) % 80 == 79)
                            cout << endl;
                        }
                    if ((x - repeat_element.first) % 80 > 0)
                        cout << endl;
                    }
                    
                    
                    element_count ++;
                }
            }
            if(VERBOSE){
                        cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"<< endl;
                        cout << "   R=" << repeat_num << "elementfreq=" << element_count << endl;
                    }
            fout1 << ">R=" << repeat_num << ", seedfreq=" << seedfreq << ", elementfreq="  << element_count << ", length=" << consensusend - consensusstart + 1 << ", Seed=";
            for (int i = 0; i < k; i++)
                fout1 << num_to_char(currentseed.second[i]);
            fout1 << endl;
            seq_type x;
            for (x = consensusstart; x <= consensusend; x++) {
                fout1 << num_to_char(consensus[x]);
                if ((x - consensusstart) % 80 == 79)
                    fout1 << endl;
                }
            if ((x - consensusstart) % 80 > 0)
                fout1 << endl;
            repeat_num++;
        }else{
            /*mask_flagのみを更新*/
            /*並列化について検討*/
            for(int i = 0; i < seedfreq; i++){
                repeat_element = masking_align(i, consensusstart, consensusend);
                maskbyrepeat_element(i, repeat_element.first, repeat_element.second, mask_flag);
            }
        }      
        
    }
    seq_type x = 0;
    if(OUTPUT_MASKED_REGION == 1){
        for (seq_type i = PADLENGTH; i < length - PADLENGTH; i++) {
            if (mask_flag2[i] == false)
                fout2 << num_to_char(sequence[i]);
            else
                fout2 << "X";
            x++;
            if (x % 50 == 0)
                fout2 << endl;
        }
    }
    freespace();
    fout1.close();
    if(OUTPUT_MASKED_REGION == 1) fout2.close();
    display_time("allrepeats");
}

void removetandem(vector<seq_type>& occs) {
    sort(occs.begin(), occs.end());
    seq_type occ;
    seq_type prevocc = -TANDEMDIST;
    auto itr = occs.begin();
    while (itr != occs.end()) {
        occ = *itr;
        if (occ < prevocc + TANDEMDIST)
            itr = occs.erase(itr);
        else
            itr++;
        prevocc = occ;
    }
}

void removemasked(vector<seq_type>& occs, const vector<bool>& mask_flag, bool isrc) {
    auto itr = occs.begin();
    int tmpflag;
    if (!isrc) {
        while (itr != occs.end()) {
            tmpflag = 0;
            for (int i = 0; i < k; i++) {
                if (mask_flag[(*itr) + i])
                    tmpflag++;
            }
            if (tmpflag == 0)
                itr++;
            else
                itr = occs.erase(itr);
        }
    } else {
        while (itr != occs.end()) {
            tmpflag = 0;
            for (int i = 0; i < k; i++) {
                if (mask_flag[(*itr) + k - 1 - i])
                    tmpflag++;
            }
            if (tmpflag == 0)
                itr++;
            else
                itr = occs.erase(itr);
        }
    }
}

void maskbyseed(const vector<seq_type>& occs, vector<bool>& mask_flag, bool isrc) {
    if (!isrc) {
        for (const auto& e : occs) {
            for (int i = 0; i < k; i++)
                mask_flag[e + i] = true;
        }
    } else {
        for (const auto& e : occs) {
            for (int i = 0; i < k; i++)
                mask_flag[e + k - 1 - i] = true;
        }
    }
}

void maskbyrepeat(int seedfreq, const vector<seq_type>& repeatstart, const vector<seq_type>& repeatend, vector<bool>& mask_flag) {
    for (int i = 0; i < seedfreq; i++) {
        if (!rev[i]) {
            for (seq_type j = pos[i] + repeatstart[i]; j <= pos[i] + repeatend[i]; j++)
                mask_flag[j] = true;
        } else {
            for (seq_type j = pos[i] - repeatend[i] + k - 1; j <= pos[i] - repeatstart[i] + k - 1; j++)
                mask_flag[j] = true;
        }
    }
}

//element startがrepeat endと対応していないことに注意(elementstartでは既に向きを考慮している)
void maskbyrepeat_element(int i, seq_type elementstart, seq_type elementend, vector<bool>& mask_flag){/*20211105_takeda_edit*/
    for (seq_type j = pos[i] + elementstart; j <= pos[i] + elementend; j++)
                mask_flag[j] = true;
}

vector<seq_type> findkmer(const vector<char>& query, const vector<vector<tuple<seq_type, seq_type, char, char>>>& cachetable) {
    int cachelength = 0;
    int tmpsz = cachetable.size();
    while (tmpsz > 3) {
        tmpsz = tmpsz >> 2;
        cachelength++;
    }
    vector<seq_type> occs;
    set<tuple<seq_type, seq_type, char, char>> matched;
    int index = 0;
    for (int i = 0; i < cachelength; i++)
        index += query[i] << i * 2;
    const auto& cache = cachetable[index];
    for (const auto& e1 : cache) {
        SA_search(query, get<0>(e1), get<1>(e1), cachelength, get<2>(e1), get<3>(e1), matched);
        for (const auto& e2 : matched) {
            for (seq_type i = get<0>(e2); i < get<1>(e2); i++)
                occs.push_back(SA[i]);
        }
        matched.clear();
    }

    return occs;
}

void SA_search(const vector<char>& query, seq_type begin, seq_type end, char query_num, char seq_num, char rem_dist, set<tuple<seq_type, seq_type, char, char>>& matched) {
    if (rem_dist >= 0 && query_num == query.size())
        matched.emplace(begin, end, seq_num, rem_dist);
    else if (rem_dist >= 0) {
        seq_type imid, imin, imax;
        seq_type b[5];
        b[0] = begin;
        for (int i = 0; i < 4; i++) {
            if (b[i] == end) {
                b[i + 1] = end;
                continue;
            }
            imin = b[i];
            imax = end - 1;
            while (1) {
                imid = (imin + imax) / 2;
                if (sequence[SA[imid] + seq_num] == i) {
                    if (imid == end - 1) {
                        b[i + 1] = end;
                        break;
                    } else if (sequence[SA[imid + 1] + seq_num] > i) {
                        b[i + 1] = imid + 1;
                        break;
                    } else
                        imin = imid + 1;
                } else {
                    if (imid == b[i]) {
                        b[i + 1] = b[i];
                        break;
                    } else
                        imax = imid - 1;
                }
            }
        }
        for (int i = 0; i < 4; i++) {
            if ((b[i + 1] - b[i]) > 0) {
                SA_search(query, b[i], b[i + 1], query_num + 1, seq_num + 1, rem_dist - (query[query_num] != i), matched);
            }
        }
    }
}

pair<int, vector<char>> find_bestseed(priority_queue<pair<int, vector<char>>>& kmers, const vector<vector<tuple<seq_type, seq_type, char, char>>>& cachetable, const vector<bool>& mask_flag) {
    int cachelength = 0;
    int tmpsz = cachetable.size();
    while (tmpsz > 3) {
        tmpsz = tmpsz >> 2;
        cachelength++;
    }

    pair<int, vector<char>> tmpbest;
    int newfreq;
    vector<char> query(k);
    vector<seq_type> occs;
    vector<seq_type> rcoccs;

    while (1) {
        tmpbest = kmers.top();
        kmers.pop();

        occs = findkmer(tmpbest.second, cachetable);
        removetandem(occs);
        removemasked(occs, mask_flag, false);
        rcoccs = findkmer(reverse_complement(tmpbest.second), cachetable);
        removetandem(rcoccs);
        removemasked(rcoccs, mask_flag, true);

        newfreq = occs.size() + rcoccs.size();

        if (kmers.empty())
            break;
        if (newfreq >= kmers.top().first)
            break;
        if (newfreq >= MINFREQ)
            kmers.push(make_pair(newfreq, tmpbest.second));
    }
    
    pos.clear();
    pos.reserve(newfreq);
    copy(occs.begin(), occs.end(), back_inserter(pos));
    copy(rcoccs.begin(), rcoccs.end(), back_inserter(pos));
    rev.resize(newfreq);
    fill(rev.begin(), rev.begin() + occs.size(), false);
    fill(rev.begin() + occs.size(), rev.end(), true);

    return make_pair(newfreq, tmpbest.second);
}

int extend(bool isright, int seedfreq, vector<seq_type>& seed_ext) {
    vector<int> nexttotalscore(4);
    vector<int> bestscore(seedfreq, 0);
    vector<vector<int>> score(seedfreq, vector<int>(2 * OFFSETWIDTH + 1, MINSCORE));//max取る用&更新
    vector<vector<int>> score_m(seedfreq, vector<int>(2 * OFFSETWIDTH + 1, MINSCORE));//match時
    vector<vector<int>> score_ins(seedfreq, vector<int>(2 * OFFSETWIDTH + 1, MINSCORE));//consensusに対するinsertion
    vector<vector<int>> score_del(seedfreq, vector<int>(2 * OFFSETWIDTH + 1, MINSCORE));//consensusに対するdeletion
    vector<vector<vector<int>>> score_m_bybase(4, vector<vector<int>>(seedfreq, vector<int>(2 * OFFSETWIDTH + 1, MINSCORE)));
    vector<vector<vector<int>>> score_ins_bybase(4, vector<vector<int>>(seedfreq, vector<int>(2 * OFFSETWIDTH + 1, MINSCORE)));
    vector<vector<vector<int>>> score_del_bybase(4, vector<vector<int>>(seedfreq, vector<int>(2 * OFFSETWIDTH + 1, MINSCORE)));
    int tmpbesttotalscore;
    int tmpbestscore;
    int bestoffset;
    char bestbase;
    int best_ext = -1;
    int besttotalscore = 0;

/*initialize*/
//20220207日報参照, (0,0)がmatchしたとすると, del方向のみスコアの初期化が必要
    for (int se = 0; se < seedfreq; se++) {
        score_m[se][OFFSETWIDTH] = 0;
        for (int offset = OFFSETWIDTH; offset > 0; offset--)
            score_del[se][offset + OFFSETWIDTH] = GAPSCORE + GAPEXTENDSCORE * (offset - 1);
    }

    for (int ext = 0; ext < MAXEXTEND; ext++) {
        fill(nexttotalscore.begin(), nexttotalscore.end(), 0);
        for (char base = 0; base < 4; base++) {
            for (int se = 0; se < seedfreq; se++)
                nexttotalscore[base] += max(0, max(bestscore[se] + CAPPENALTY, compute_score(isright, ext, se, base, score_m, score_ins, score_del, score_m_bybase, score_ins_bybase, score_del_bybase)));
        }
        bestbase = distance(nexttotalscore.begin(), max_element(nexttotalscore.begin(), nexttotalscore.end()));

        if (isright)
            consensus[MAXEXTEND + ext] = bestbase;
        else
            consensus[MAXEXTEND - ext - 1] = bestbase;

        for (int se = 0; se < seedfreq; se++) {
            for (int offset = -OFFSETWIDTH; offset <= OFFSETWIDTH; offset++) {
                score_m[se][offset + OFFSETWIDTH] = score_m_bybase[bestbase][se][offset + OFFSETWIDTH];
                score_ins[se][offset + OFFSETWIDTH] = score_ins_bybase[bestbase][se][offset + OFFSETWIDTH];
                score_del[se][offset + OFFSETWIDTH] = score_del_bybase[bestbase][se][offset + OFFSETWIDTH];
                score[se][offset + OFFSETWIDTH] = max(score_m[se][offset + OFFSETWIDTH], max(score_ins[se][offset + OFFSETWIDTH], score_del[se][offset + OFFSETWIDTH]));
            }
        }
        tmpbesttotalscore = 0;
        for (int se = 0; se < seedfreq; se++) {
            tmpbestscore = MINSCORE * 2;
            for (int offset = -OFFSETWIDTH; offset <= OFFSETWIDTH; offset++) {
                if (score[se][offset + OFFSETWIDTH] > tmpbestscore) {
                    tmpbestscore = score[se][offset + OFFSETWIDTH];
                    bestoffset = offset;
                }
            }
            if (tmpbestscore > bestscore[se]) {
                seed_ext[se] = ext + bestoffset;
                bestscore[se] = tmpbestscore;
            }
            tmpbesttotalscore += max(tmpbestscore, bestscore[se] + CAPPENALTY);
        }

        if (tmpbesttotalscore >= besttotalscore + (ext - best_ext) * MINIMPROVEMENT) {
            best_ext = ext;
            besttotalscore = tmpbesttotalscore;
        }

        //cout << "ext:" << ext << "\tbest_ext:" << best_ext << endl;
        if (ext - best_ext >= WHEN_TO_STOP)
            break;
    }

    return best_ext;
}

int compute_score(bool isright, int ext, int se, char base, const vector<vector<int>>& score_m, const vector<vector<int>>& score_ins, const vector<vector<int>>& score_del, vector<vector<vector<int>>>& score_m_bybase, vector<vector<vector<int>>>& score_ins_bybase, vector<vector<vector<int>>>& score_del_bybase) {
    vector<int> nextscore_m_bybase(2 * OFFSETWIDTH + 1);
    vector<int> nextscore_ins_bybase(2 * OFFSETWIDTH + 1,MINSCORE);
    vector<int> nextscore_del_bybase(2 * OFFSETWIDTH + 1,MINSCORE);
    int tmpscore = MINSCORE;

    /*nextscore_m_bybase*/
    if (!rev[se] && isright) {
        for (int offset = -OFFSETWIDTH; offset <= OFFSETWIDTH; offset++) {
            if (base == sequence[pos[se] + offset + ext])
                nextscore_m_bybase[offset + OFFSETWIDTH] = max(score_m[se][offset + OFFSETWIDTH], max(score_ins[se][offset + OFFSETWIDTH], score_del[se][offset + OFFSETWIDTH])) + MATCHSCORE;
            else
                nextscore_m_bybase[offset + OFFSETWIDTH] = max(score_m[se][offset + OFFSETWIDTH], max(score_ins[se][offset + OFFSETWIDTH], score_del[se][offset + OFFSETWIDTH])) + MISMATCHSCORE;
        }
    } else if (rev[se] && isright) {
        for (int offset = -OFFSETWIDTH; offset <= OFFSETWIDTH; offset++) {
            if (base == complement(sequence[pos[se] + k - 1 - (offset + ext)]))
                nextscore_m_bybase[offset + OFFSETWIDTH] = max(score_m[se][offset + OFFSETWIDTH], max(score_ins[se][offset + OFFSETWIDTH], score_del[se][offset + OFFSETWIDTH])) + MATCHSCORE;
            else
                nextscore_m_bybase[offset + OFFSETWIDTH] = max(score_m[se][offset + OFFSETWIDTH], max(score_ins[se][offset + OFFSETWIDTH], score_del[se][offset + OFFSETWIDTH])) + MISMATCHSCORE;
        }
    } else if (!rev[se] && !isright) {
        for (int offset = -OFFSETWIDTH; offset <= OFFSETWIDTH; offset++) {
            if (base == sequence[pos[se] - 1 - offset - ext])
                nextscore_m_bybase[offset + OFFSETWIDTH] = max(score_m[se][offset + OFFSETWIDTH], max(score_ins[se][offset + OFFSETWIDTH], score_del[se][offset + OFFSETWIDTH])) + MATCHSCORE;
            else
                nextscore_m_bybase[offset + OFFSETWIDTH] = max(score_m[se][offset + OFFSETWIDTH], max(score_ins[se][offset + OFFSETWIDTH], score_del[se][offset + OFFSETWIDTH])) + MISMATCHSCORE;
        }
    } else {
        for (int offset = -OFFSETWIDTH; offset <= OFFSETWIDTH; offset++) {
            if (base == complement(sequence[pos[se] + k - 1 - (-1 - offset - ext)]))
                nextscore_m_bybase[offset + OFFSETWIDTH] = max(score_m[se][offset + OFFSETWIDTH], max(score_ins[se][offset + OFFSETWIDTH], score_del[se][offset + OFFSETWIDTH])) + MATCHSCORE;
            else
                nextscore_m_bybase[offset + OFFSETWIDTH] = max(score_m[se][offset + OFFSETWIDTH], max(score_ins[se][offset + OFFSETWIDTH], score_del[se][offset + OFFSETWIDTH])) + MISMATCHSCORE;
        }
    }

    /*nextscore_ins_bybase*/
    //score_ins_bybase[offset + OFFSETWIDTH] = MINSCORE;
    for (int offset = -OFFSETWIDTH; offset < OFFSETWIDTH; offset++) {
        nextscore_ins_bybase[offset + OFFSETWIDTH] = max(score_m[se][offset + 1 + OFFSETWIDTH] + GAPSCORE, score_ins[se][offset + 1 + OFFSETWIDTH] + GAPEXTENDSCORE);
    }

    /*nextscore_del_bybase*/
    //score_del_bybase[-offset + OFFSETWIDTH] = MINSCORE;
    for (int offset = -OFFSETWIDTH + 1; offset <= OFFSETWIDTH; offset++) {
        nextscore_del_bybase[offset + OFFSETWIDTH] = max(nextscore_m_bybase[offset -1 + OFFSETWIDTH] + GAPSCORE, nextscore_del_bybase[offset -1 + OFFSETWIDTH] + GAPEXTENDSCORE);
    }

    /*2022_02_07 最大値取るところから,スコア用の行列用意しても良いかも*/

    for (int offset = -OFFSETWIDTH; offset <= OFFSETWIDTH; offset++) {
        score_m_bybase[base][se][offset + OFFSETWIDTH] = nextscore_m_bybase[offset + OFFSETWIDTH];
        score_ins_bybase[base][se][offset + OFFSETWIDTH] = nextscore_ins_bybase[offset + OFFSETWIDTH];
        score_del_bybase[base][se][offset + OFFSETWIDTH] = nextscore_del_bybase[offset + OFFSETWIDTH];
        tmpscore = max(tmpscore,max(nextscore_m_bybase[offset + OFFSETWIDTH], max(nextscore_ins_bybase[offset + OFFSETWIDTH], nextscore_del_bybase[offset + OFFSETWIDTH])));
    }

    return tmpscore;
}

pair<seq_type, seq_type> masking_align(int i,seq_type consensusstart, seq_type consensusend){
    vector<int> right_mask_score(2 * OFFSETWIDTH + 1, MINSCORE);
    vector<int> left_mask_score(2 * OFFSETWIDTH + 1, MINSCORE);
    vector<int> best_right_mask_score(2 * OFFSETWIDTH + 1, MINSCORE);
    vector<int> best_left_mask_score(2 * OFFSETWIDTH + 1, MINSCORE);

    vector<int> right_mask_score_m(2 * OFFSETWIDTH + 1, MINSCORE);//match時
    vector<int> right_mask_score_ins(2 * OFFSETWIDTH + 1, MINSCORE);//consensusに対するinsertion
    vector<int> right_mask_score_del(2 * OFFSETWIDTH + 1, MINSCORE);//consensusに対するdeletion

    vector<int> left_mask_score_m(2 * OFFSETWIDTH + 1, MINSCORE);//match時
    vector<int> left_mask_score_ins(2 * OFFSETWIDTH + 1, MINSCORE);//consensusに対するinsertion
    vector<int> left_mask_score_del(2 * OFFSETWIDTH + 1, MINSCORE);//consensusに対するdeletion
    
    int right_ext, left_ext, right_best_ext = MAXEXTEND, left_best_ext = MAXEXTEND - 1;

    int right_tmp_score,left_tmp_score;
    int right_bestscore = 0, left_bestscore = 0;

    int elementstart, elementend;

    /*right extend*/
    right_mask_score_m[OFFSETWIDTH] = 0;
    for (int offset = OFFSETWIDTH; offset > 0; offset--)
        right_mask_score_del[offset + OFFSETWIDTH] = GAPSCORE + GAPEXTENDSCORE * (offset - 1);

    for(right_ext = MAXEXTEND; right_ext <= consensusend; right_ext ++){
        right_tmp_score = mask_extention_score(true, right_ext, i, right_mask_score, right_mask_score_m, right_mask_score_ins,  right_mask_score_del);//cap penarty入れる必要あるかも, 要修正
        if(right_tmp_score > right_bestscore){
            right_bestscore = right_tmp_score;
            right_best_ext = right_ext;
            best_right_mask_score = right_mask_score;
        }

        if(right_ext -right_best_ext >= WHEN_TO_STOP) break;
    }
    /*left extend*/
    left_mask_score_m[OFFSETWIDTH] = 0;
    for (int offset = OFFSETWIDTH; offset > 0; offset--)
        left_mask_score_del[offset + OFFSETWIDTH] = GAPSCORE + GAPEXTENDSCORE * (offset - 1);

    for(left_ext = MAXEXTEND - 1; left_ext >= consensusstart; left_ext --){
        left_tmp_score = mask_extention_score(false, left_ext, i, left_mask_score, left_mask_score_m, left_mask_score_ins, left_mask_score_del);//cap penarty入れる必要あるかも, 要修正
        if(left_tmp_score > left_bestscore){
            left_bestscore = left_tmp_score;
            left_best_ext = left_ext;
            best_left_mask_score = left_mask_score;
        }
        //if(VERBOSE) cout << left_tmp_score <<endl;
        if(left_best_ext - left_ext >= WHEN_TO_STOP) break;
    }
    /*curation of extension result*/
    int best_rightoffset = -OFFSETWIDTH;
    int best_leftoffset = -OFFSETWIDTH;

    for(int rightoffset = -OFFSETWIDTH; rightoffset <= OFFSETWIDTH; rightoffset ++){
        if(best_right_mask_score[best_rightoffset + OFFSETWIDTH ] < best_right_mask_score[rightoffset + OFFSETWIDTH]){
            best_rightoffset = rightoffset;
        }
    }
    for(int leftoffset = -OFFSETWIDTH; leftoffset <= OFFSETWIDTH; leftoffset ++){
        if(best_left_mask_score[best_leftoffset + OFFSETWIDTH] < best_left_mask_score[leftoffset + OFFSETWIDTH]){
            best_leftoffset = leftoffset;
        }
    }
   
   if(VERBOSE == true){
       printf("best_leftoffset=%d\n",best_leftoffset);
       printf("best_rightoffset=%d\n",best_rightoffset);
   }
   if(!rev[i]){
       elementstart = left_best_ext - MAXEXTEND - best_leftoffset;
       elementend = right_best_ext - MAXEXTEND + best_rightoffset;
   }else{
       elementstart = - right_best_ext + k - 1 + MAXEXTEND - best_rightoffset;
       elementend =  - left_best_ext + k - 1  + MAXEXTEND + best_leftoffset;//20211222 revのleft_extが1塩基分多く数えていて, さらに"-1"する必要がある可能性アリ, 要検討
   }
    printf("elementstart=%d\n",elementstart);
    printf("elementend=%d\n",elementend);
    return make_pair(elementstart, elementend);
}


int mask_extention_score(bool isright, int ext, int i, vector<int> &mask_score, vector<int> &mask_score_m, vector<int> &mask_score_ins, vector<int> &mask_score_del){//20211222 デバッグ完了, 合ってそう, 20220208 affine gap
    vector<int> mask_nextscore_m(2 * OFFSETWIDTH + 1);
    vector<int> mask_nextscore_ins(2 * OFFSETWIDTH + 1);
    vector<int> mask_nextscore_del(2 * OFFSETWIDTH + 1);
    
    int mask_tmpscore = MINSCORE;

    if (!rev[i] && isright) {
        for (int offset = -OFFSETWIDTH; offset <= OFFSETWIDTH; offset++) {
            if (consensus[ext] == sequence[pos[i] + offset + ext - MAXEXTEND]){
                mask_nextscore_m[offset + OFFSETWIDTH] = max(mask_score_m[offset + OFFSETWIDTH], max(mask_score_ins[offset + OFFSETWIDTH], mask_score_del[offset + OFFSETWIDTH])) + MATCHSCORE;
            }
            else{
                mask_nextscore_m[offset + OFFSETWIDTH] = max(mask_score_m[offset + OFFSETWIDTH], max(mask_score_ins[offset + OFFSETWIDTH], mask_score_del[offset + OFFSETWIDTH])) + MISMATCHSCORE;
            }
        }     

    }else if (rev[i] && isright) {
        for (int offset = -OFFSETWIDTH; offset <= OFFSETWIDTH; offset++) {
            if (consensus[ext] == complement(sequence[pos[i] + k - 1 - (- offset + ext - MAXEXTEND)])){
                mask_nextscore_m[offset + OFFSETWIDTH] = max(mask_score_m[offset + OFFSETWIDTH], max(mask_score_ins[offset + OFFSETWIDTH], mask_score_del[offset + OFFSETWIDTH])) + MATCHSCORE;
            }
            else{
                mask_nextscore_m[offset + OFFSETWIDTH] = max(mask_score_m[offset + OFFSETWIDTH], max(mask_score_ins[offset + OFFSETWIDTH], mask_score_del[offset + OFFSETWIDTH])) + MISMATCHSCORE;
            }
        }
        
    }else if (!rev[i] && !isright) {
        for (int offset = -OFFSETWIDTH; offset <= OFFSETWIDTH; offset++) {
            if (consensus[ext] == sequence[pos[i] - offset + (ext - MAXEXTEND)]){
                mask_nextscore_m[offset + OFFSETWIDTH] = max(mask_score_m[offset + OFFSETWIDTH], max(mask_score_ins[offset + OFFSETWIDTH], mask_score_del[offset + OFFSETWIDTH])) + MATCHSCORE;
            }
            else{
                mask_nextscore_m[offset + OFFSETWIDTH] = max(mask_score_m[offset + OFFSETWIDTH], max(mask_score_ins[offset + OFFSETWIDTH], mask_score_del[offset + OFFSETWIDTH])) + MISMATCHSCORE;
            }
        }

    }else{
        for (int offset = -OFFSETWIDTH; offset <= OFFSETWIDTH; offset++) {
            if (consensus[ext] == complement(sequence[pos[i] + k - 1 - (offset + ext - MAXEXTEND)])){
                mask_nextscore_m[offset + OFFSETWIDTH] = max(mask_score_m[offset + OFFSETWIDTH], max(mask_score_ins[offset + OFFSETWIDTH], mask_score_del[offset + OFFSETWIDTH])) + MATCHSCORE;
            }
            else{
                mask_nextscore_m[offset + OFFSETWIDTH] = max(mask_score_m[offset + OFFSETWIDTH], max(mask_score_ins[offset + OFFSETWIDTH], mask_score_del[offset + OFFSETWIDTH])) + MISMATCHSCORE;
            }
        }
    }

    for (int offset = -OFFSETWIDTH; offset < OFFSETWIDTH; offset++) {
        mask_nextscore_ins[offset + OFFSETWIDTH] = max(mask_score_m[offset + 1 + OFFSETWIDTH] + GAPSCORE, mask_score_ins[offset + 1 + OFFSETWIDTH] + GAPEXTENDSCORE);
    }

    for (int offset = -OFFSETWIDTH + 1; offset <= OFFSETWIDTH; offset++) {
        mask_nextscore_del[offset + OFFSETWIDTH] = max(mask_nextscore_m[offset -1 + OFFSETWIDTH] + GAPSCORE, mask_nextscore_del[offset -1 + OFFSETWIDTH] + GAPEXTENDSCORE);
    }

    for (int offset = -OFFSETWIDTH; offset <= OFFSETWIDTH; offset++) {
        mask_score_m[offset + OFFSETWIDTH] = mask_nextscore_m[offset + OFFSETWIDTH];
        mask_score_ins[offset + OFFSETWIDTH] = mask_nextscore_ins[offset + OFFSETWIDTH];
        mask_score_del[offset + OFFSETWIDTH] = mask_nextscore_del[offset + OFFSETWIDTH];
        mask_score[offset + OFFSETWIDTH] = max(mask_score_m[offset + OFFSETWIDTH], max(mask_score_ins[offset + OFFSETWIDTH], mask_score_del[offset + OFFSETWIDTH]));
        mask_tmpscore = max(mask_tmpscore,mask_score[offset + OFFSETWIDTH]);
    }

    return mask_tmpscore;
}

/*2021105_chrtableについて編集*/
/*20220207_chrtableにpaddeingを入れる処置をとった*/
void build_sequence() {
    ifstream fin;
    char c;
    string str;
    seq_type MAXLENGTH;
    int NUMCHR = 0;

    fin.open(SEQUENCE_FILE);
    if (!fin) {
        cerr << "Could not open sequence file " << SEQUENCE_FILE << endl;
        exit(1);
    }
    fin.seekg(0, ios_base::end);
    MAXLENGTH = fin.tellg();
    fin.seekg(0, ios_base::beg);


    while (getline(fin,str)){
        if (str[0] == '>') NUMCHR ++;
    }
    fin.clear();
    fin.seekg(0, ios_base::beg);
    //debug
    cout << "NUMCHR " << NUMCHR << endl;

    sequence = (char*)malloc((MAXLENGTH + PADLENGTH * (2 + NUMCHR)) * sizeof(char));
    if (NULL == sequence) {
        cerr << "Unable to allocate " << MAXLENGTH + k * 2 << " bytes for the sequence" << endl;
        exit(1);
    }

    length = 0;
    chrtable.push_back(make_pair("unknown",0));

    for (int i = 0; i < PADLENGTH; i++) {
        sequence[length] = 99;
        length++;
    }
    while (!fin.eof()) {
        c = fin.get();

        if (c == EOF)
            break;
        if (c == '\n')
            continue;
        if (c == '>') {
            chrtable.push_back(make_pair("padding",length));
            for (int i = 0; i < PADLENGTH; i++) {
                sequence[length] = 99;
                length++;
            }     
            
            string chrname = "";
            while (!fin.eof()) {
                fin.get(c);
                if (c == '\n'){
                    chrtable.push_back(make_pair(chrname,length));
                    break;
                }
                if (c == '\t' || c == ' '){
                    chrname += '_';
                }else{
                chrname += c;
                }
            }
        } else {
            while (!fin.eof()) {
                if (c > 64) {
                    sequence[length] = char_to_num(c);
                    length++;
                }
                fin.get(c);
                if (c == '\n')
                    break;
            }
        }
    }
    for (int i = 0; i < PADLENGTH; i++) {
        sequence[length] = 99;
        length++;
    }
    fin.close();
}

pair<string, seq_type> chrtracer(const seq_type stringpos){
    pair<string, seq_type> output = chrtable[int(chrtable.size()-1)];
    for(int i = 0; i < chrtable.size() ; i++){
        if(stringpos < chrtable[i].second){
            output = chrtable[i-1];
            break;
        }     
    }
    return output;
}

void allocate_space() {
    consensus = (char*)malloc((2 * MAXEXTEND) * sizeof(char));
    consensus[2 * MAXEXTEND] = (char)NULL;
    if (NULL == consensus) {
        cerr << "Could not allocate space for consensus array" << endl;
        exit(1);
    }
}

void freespace() {
    free(consensus);
    //free(rev);
    //free(pos);
}

char num_to_char(char z) {
    if (z == 0)
        return (char)'A';
    if (z == 1)
        return (char)'C';
    if (z == 2)
        return (char)'G';
    if (z == 3)
        return (char)'T';
    return (char)'N';
}

char char_to_num(char c) {
    if (c == 'A' || c == 'a')
        return 0;
    if (c == 'C' || c == 'c')
        return 1;
    if (c == 'G' || c == 'g')
        return 2;
    if (c == 'T' || c == 't')
        return 3;
    if (c == 'N' || c == 'n')
        return 99;
    if (c == 'x')
        return 99;
    if (IUPAC(c))
        return 99;
    cerr << "char_to_num error: c=" << c << endl;
    exit(1);
}

char complement(char c) {
    if (c == 0)
        return 3;
    if (c == 1)
        return 2;
    if (c == 2)
        return 1;
    if (c == 3)
        return 0;
    return 99;
}

vector<char> reverse_complement(const vector<char>& query) {
    vector<char> ret(k);
    for (int i = 0; i < k; i++)
        ret[i] = complement(query[k - i - 1]);
    return ret;
}

double compute_entropy(const vector<char>& kmer) {
    int x;
    int count[4];
    double answer, y;

    for (x = 0; x < 4; x++)
        count[x] = 0;
    for (x = 0; x < k; x++)
        count[(int)kmer[x]] += 1;
    answer = 0.0;
    for (x = 0; x < 4; x++) {
        if (count[x] == 0)
            continue;
        y = ((double)count[x]) / ((double)k);
        answer += y * log(y);
    }
    return answer;
}

/*
 * Set k = ceil( 1 + log(\sum(g.segments[*].length)/log(4) );
 */
int default_k(seq_type len, int KMERDIST) { //20210405変更
    vector<vector<long long int> > v(40,vector<long long int>(40));//パスカルの三角形
    for(int i = 0;i <v.size(); i++){
        v[i][0]=1;
        v[i][i]=1;
    }
    for(int kk = 1;kk <v.size();kk++){ //kを使っていたのでkkにしただけ. 
        for(int j = 1;j<kk;j++){
            v[kk][j]=(v[kk-1][j-1]+v[kk-1][j]);
        }
    }

    int l = 2;
    double E;
    for(; l <= 40; l ++){
        double Comb = 0;
        for(int d = 0; d <= KMERDIST; d ++){
            Comb += v[l][d] * pow(3,d);
        }
        E = len * Comb/pow(4,l);
        if(E<1) break;
    }

    return l + 1;
}

void display_time(string msg) {
#ifdef _OPENMP
    double time0 = time1;
    time1 = omp_get_wtime();
    totaltime += time1 - time0;
    cout << msg << ":" << round(time1 - time0) << "sec Total:" << round(totaltime) << "sec";
#else
    clock_t time0 = time1;
    time1 = clock();
    totaltime += time1 - time0;
    cout << msg << ":" << round((double)(time1 - time0) / CLOCKS_PER_SEC) << "sec Total:" << round((double)totaltime / CLOCKS_PER_SEC) << "sec";
#endif
    struct rusage r;
    getrusage(RUSAGE_SELF, &r);
    cout << " maxrss=" << r.ru_maxrss << "KB" << endl;
}


void print_usage(){
    cout << "REPrise: de novo interspersed repeat detection software. version 1.0.1" << endl;
    cout << endl;
    cout << "Usage" << endl;
    cout << endl;
    cout << "REPrise [-input genome file] [-output outputname] [Options]" << endl;
    cout << endl;
    cout << "Options\n";
    cout << "(Required)\n";
    cout << "   -input  STR         input file name. You can input assembled genome file, or hard masked genome file\n";
    cout << "   -output STR         output file name. REPrise outputs STR.freq, STR.bed STR.masked and STR.reprof (consensus seqnences)\n" ;
    cout << endl;
    cout << "(Optional)\n";
    cout << "   -h                  Print help and exit\n";
    cout << "   -v                  Verbose\n";
    cout << "   -additonalfile      Output files about masked region(.masked and .bed)\n";
    cout << endl;
    cout << "   -match INT          Match score of the extension alignment (default = 1)\n";
    cout << "   -match INT          Mismatch score of the extension alignment (default = -1)\n";
    cout << "   -gap   INT          Gap open score of the extension alignment (default = -5)\n";
    cout << "   -gapex  INT         Gap extension score of the extension alignment (default = -1)\n";
    cout << "   -capplenalty INT    Penalty of the imcomplete length alignment (default = -20)\n";
    cout << "   -dist INT           Number of mismatches allowed in inexact seed (default = 0)\n";
    cout << endl;
    cout << "   -maxextend INT      Uppler limit length of extension in one side direction of consensus repeat (default = 10000)\n";
    cout << "   -maxrepeat INT      Maximum Number of elements belonging to one repeat family (default = 100000)\n";
    cout << "   -maxgap INT         Band size(= maximum number of gaps allowed) of extension alignment (default = 5)\n";
    cout << "   -stopafter INT      If the maximum score of extension alignment does not change INT consecutive times, that alignment will stop (default = 100)\n";
    cout << "   -minlength INT      Minimum number of length of the consensus sequence of repeat family(default = 50)\n";
    cout << "   -minfreq INT        Minimum number of elements  belonging to one repeat family (default = 3)\n";
    cout << "   -minimprovement INT Penalty associated with the number of regions to be extended as the repeat regions (default = 3)\n";
    cout << "   -tandemdist INT     Interval to match the same seed to avoid seed matching with tandem repeats(default = 500)\n";
    cout << endl;
    cout << "   -pa INT             Number of openMP parallel cores";

    cout << endl;
}