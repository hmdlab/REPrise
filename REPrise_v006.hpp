#ifndef __REPrise_H__
#define __REPrise_H__

#define seq_type long

void store_cache(int edit_distance, vector<vector<tuple<seq_type, seq_type, char, char>>> &cachetable);
void build_sortedkmers(priority_queue<pair<int, vector<char>>> &kmers, const vector<vector<tuple<seq_type, seq_type, char, char>>> &dist0cachetable, const vector<vector<tuple<seq_type, seq_type, char, char>>> &cachetable);
void build_repeat_families(priority_queue<pair<int, vector<char>>> &kmers, const vector<vector<tuple<seq_type, seq_type, char, char>>> &cachetable);

vector<seq_type> findkmer(const vector<char>& query, const vector<vector<tuple<seq_type, seq_type, char, char>>> &cachetable);
void SA_search(const vector<char> &query, seq_type begin, seq_type end, char query_num, char seq_num, char rem_dist, set<tuple<seq_type, seq_type, char, char>> &matched);
pair<int, vector<char>> find_bestseed(priority_queue<pair<int, vector<char>>> &kmers, const vector<vector<tuple<seq_type, seq_type, char, char>>> &cachetable, const vector<bool> &mask_flag);
int extend(bool isright, int seedfreq, vector<seq_type> &seed_ext);
int compute_score(bool isright, int ext, int se, char base, const vector<vector<int>>& score_m, const vector<vector<int>>& score_ins, const vector<vector<int>>& score_del, vector<vector<vector<int>>>& score_m_bybase, vector<vector<vector<int>>>& score_ins_bybase, vector<vector<vector<int>>>& score_del_bybase);

void removetandem(vector<seq_type> &occs);
void removemasked(vector<seq_type> &occs, const vector<bool> &mask_flag, bool isrc);
void maskbyseed(const vector<seq_type> &occs, vector<bool> &mask_flag, bool isrc);
void maskbyrepeat(int seedfreq, const vector<seq_type> &repeatstart, const vector<seq_type> &repeatend, vector<bool> &mask_flag);
void maskbyrepeat_element(int i, seq_type elementstart, seq_type elementend, vector<bool> &mask_flag);/*20211105_takeda_edit*/

/*20211105_takeda_edit*/
pair<seq_type, seq_type> masking_align(int i,seq_type consensusstart, seq_type consensusend);
int mask_extention_score(bool isright, int ext, int i, vector<int> &mask_score, vector<int> &mask_score_m, vector<int> &mask_score_ins, vector<int> &mask_score_del);

void build_sequence();
pair<string, seq_type> chrtracer(const seq_type stringpos);
void allocate_space();
void freespace();
char num_to_char(char z);
char char_to_num(char c);
char complement(char c);
vector<char> reverse_complement(const vector<char> &query);
double compute_entropy(const vector<char> &kmer);
int default_k(seq_type len,int KMERDIST);

void display_time(string msg);
#endif
