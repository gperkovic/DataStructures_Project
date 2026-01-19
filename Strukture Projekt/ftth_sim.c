/*
  Naziv programa:  FTTH/GPON Simulator

Program:
    - čita datoteku sa topologijom optičke mreže (OLT-SPLITTER-ONT)
    - računa optičke gubitke i RX snagu za svaki ONT (za svakog korisnika)
    - sakuplja i grupira statistiku po splitteru
    - stvara datoteke ont_results.csv i splitter_results.csv iz kojih dobivamo vizualini prikaz pomoću grafova
    - stvara datoteku report.txt kojom izvještava status optičke mreže
*/

#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#define MAX_ONT 1024
#define TOP_N   5

typedef enum { NODE_OLT, 
    NODE_SPLITTER, 
    NODE_ONT 
} NodeType;

typedef struct Node {
    NodeType type;

    // Link parametri od roditelja
    double len_km;
    int connectors;
    int splices;

    // Splitter parametri
    int splitter_ratio;         // npr. 8/16/32/64 dijelitelja
    char name[64];              

    // ONT parametri
    int ont_id;

    // OLT parametri
    double olt_tx_dbm;          // snaga pokretanja
    double gpon_rxmin_dbm;      // minimalna RX snaga (prag)

    // Greška prijenosa
    int faulty;                 // ako je 1, dodaje extra_loss_db i označava granu "down-like"
    double extra_loss_db;       // dodatno gubljenje optičkog signala

    struct Node* child;
    struct Node* sibling;
} Node;

typedef struct {
    int ont_count;
    int ok_count;
    int fail_count;  // ispod RX min
    int down_count;  // prisiljen na down zbog krivog patha
    double sum_rx;
    double sum_loss;
    double best_rx;
    double worst_rx;
} SubtreeStats;

typedef struct {
    char name[64];
    int ratio;
    int ont_count;
    int ok_count;
    int fail_count;
    int down_count;
    double avg_rx;
    double avg_loss;
    double worst_rx;
} SplitterRecord;

typedef struct {
    SplitterRecord* arr;
    size_t n;
    size_t cap;
} SplitterList;

typedef struct {
    int ont_id;
    double rx_dbm;
    double margin_db;
    char path[512];
} OntResult;

OntResult ont_results[MAX_ONT];
int ont_results_count = 0;

// --------- constante optičke mreže ----------
static const double ATTEN_DB_PER_KM = 0.35;   
static const double CONN_LOSS_DB    = 0.50;
static const double SPLICE_LOSS_DB  = 0.10;
static const double SPLITTER_INS_DB = 1.00;   // insertion loss (uz idealni 10log10)
// ------------------------------------------------

static void die(const char* msg);
static Node* node_new(NodeType t);
static void node_add_child(Node* parent, Node* child);
static void splitter_list_init(SplitterList* sl);
static void splitter_list_push(SplitterList* sl, const SplitterRecord* rec);
static char* ltrim(char* s);
static void rtrim(char* s);
static int leading_spaces(const char* s);
static NodeType parse_type(const char* tok);
static int parse_int(const char* v);
static double parse_double(const char* v);
static void apply_kv(Node* n, const char* key, const char* val);
static void parse_line(Node* n, char* line_no_indent);
static double splitter_loss_db(int ratio);
static double node_link_loss_db(const Node* n);
static SubtreeStats stats_init(void);
static void stats_merge(SubtreeStats* a, const SubtreeStats* b);
static void path_append(char* path, size_t cap, const char* part);
static SubtreeStats walk_and_compute(const Node* n, double parent_tx_dbm, double rxmin_dbm, double acc_loss_db,
    double acc_dist_km, int down_flag, FILE* ont_csv, SplitterList* splitters, char* path, size_t path_cap);
static void write_splitter_csv(const char* filename, const SplitterList* sl);
static void print_summary(const SubtreeStats* all, double tx, double rxmin);
int cmp_margin(const void* a, const void* b);
static Node* read_topology(const char* filename);
void generate_report(const SubtreeStats* stats, Node* root);
static void free_tree(Node* n);

int main(int argc, char** argv) {
    if (argc < 2) {
        printf("Koristimo %s ftth_topology.txt\n", argv[0]);
        return 1;
    }

    const char* topo_file = argv[1];
    Node* root = read_topology(topo_file);

    FILE* ont_csv = fopen("ont_results.csv", "w");
    if (!ont_csv) {
        die("Nemoguce je otvoriti ont_results.csv za pisanje");
    }
    fprintf(ont_csv, "ont_id,total_dist_km,total_loss_db,rx_dbm,margin_db,status,path\n");

    SplitterList splitters;
    splitter_list_init(&splitters);

    char path[512] = {0};

    // pocinjemo rekurziju od root
    SubtreeStats all = walk_and_compute(
        root,
        root->olt_tx_dbm,
        root->gpon_rxmin_dbm,
        0.0, 0.0, 0,
        ont_csv,
        &splitters,
        path, sizeof(path)
    );

    fclose(ont_csv);

    write_splitter_csv("splitter_results.csv", &splitters);

    print_summary(&all, root->olt_tx_dbm, root->gpon_rxmin_dbm);

    qsort(ont_results, ont_results_count, sizeof(OntResult), cmp_margin);

    printf("\nTOP %d najgorih ONT-ova (po margin):\n", TOP_N);
    for (int i = 0; i < TOP_N && i < ont_results_count; i++) {
        printf(
            "ONT %d | margin = %.2f dB | RX = %.2f dBm | path = %s\n",
            ont_results[i].ont_id,
            ont_results[i].margin_db,
            ont_results[i].rx_dbm,
            ont_results[i].path
        );
    }

    printf("\nStvorene datoteke:\n");
    printf(" - ont_results.csv\n");
    printf(" - splitter_results.csv\n");

    generate_report(&all, root);
    printf("\nStvoren report.txt\n");

    free(splitters.arr);
    free_tree(root);
    return 0;
}

static void die(const char* msg) {
    fprintf(stderr, "Error: %s\n", msg);
    exit(1);
}

static Node* node_new(NodeType t) {
    Node* n = (Node*)calloc(1, sizeof(Node));
    if (!n) {
        die("Nema slobodne memorije");
    }
    n->type = t;
    n->splitter_ratio = 0;
    n->ont_id = -1;
    n->olt_tx_dbm = 3.0;
    n->gpon_rxmin_dbm = -27.0;
    n->faulty = 0;
    n->extra_loss_db = 0.0;
    n->name[0] = '\0';
    return n;
}

static void node_add_child(Node* parent, Node* child) {
    if (!parent->child) {
        parent->child = child;
        return;
    }
    Node* cur = parent->child;
    while (cur->sibling) {
        cur = cur->sibling;
    }
    cur->sibling = child;
}

static void splitter_list_init(SplitterList* sl) {
    sl->arr = NULL; 
    sl->n = 0; 
    sl->cap = 0;
}

static void splitter_list_push(SplitterList* sl, const SplitterRecord* rec) {
    if (sl->n == sl->cap) {
        size_t newcap = sl->cap ? sl->cap * 2 : 16;
        SplitterRecord* p = (SplitterRecord*)realloc(sl->arr, newcap * sizeof(SplitterRecord));

        if (!p) {
            die("Nema slobodne memorije (realloc)");
        }
        sl->arr = p;
        sl->cap = newcap;
    }
    sl->arr[sl->n++] = *rec;
}

// brisanje razmaka s lijeve strane stringa
static char* ltrim(char* s) {
    while (*s && isspace((unsigned char)*s)) {
        s++;
    }
    return s;
}

// brisanje razmaka s desne strane stringa
static void rtrim(char* s) {
    size_t n = strlen(s);
    while (n && isspace((unsigned char)s[n - 1])) {
        s[--n] = '\0';
    }
}

// broji razmake (određuje dubinu stabla)
static int leading_spaces(const char* s) {
    int c = 0;
    while (*s == ' ') { 
        c++; 
        s++; 
    }
    return c;
}

// uvlaka mora biti paran broj razmaka (2 po razini - određuje dubinu stabla)
static NodeType parse_type(const char* tok) {
    if (strcmp(tok, "OLT") == 0) {
        return NODE_OLT;
    }
    if (strcmp(tok, "SPLITTER") == 0) {
        return NODE_SPLITTER;
    }
    if (strcmp(tok, "ONT") == 0) {
        return NODE_ONT;
    }
    die("Nepoznati node type (ocekivani su OLT/SPLITTER/ONT)");
    return NODE_ONT;
}

static int parse_int(const char* v) { 
    return (int)strtol(v, NULL, 10); 
}

static double parse_double(const char* v) { 
    return strtod(v, NULL); 
}

// key=value
static void apply_kv(Node* n, const char* key, const char* val) {
    if (strcmp(key, "len") == 0) {
        n->len_km = parse_double(val);
    } else if (strcmp(key, "km") == 0) {
        n->len_km = parse_double(val);
    } else if (strcmp(key, "conn") == 0) {
        n->connectors = parse_int(val);
    } else if (strcmp(key, "sp") == 0) {
        n->splices = parse_int(val);
    } else if (strcmp(key, "ratio") == 0) {
        n->splitter_ratio = parse_int(val);
    } else if (strcmp(key, "name") == 0) {
        strncpy(n->name, val, sizeof(n->name) - 1);
        n->name[sizeof(n->name) - 1] = '\0';
    } else if (strcmp(key, "id") == 0) {
        n->ont_id = parse_int(val);
    } else if (strcmp(key, "tx") == 0) {
        n->olt_tx_dbm = parse_double(val);
    } else if (strcmp(key, "rxmin") == 0) {
        n->gpon_rxmin_dbm = parse_double(val);
    } else if (strcmp(key, "faulty") == 0) {
        n->faulty = parse_int(val) ? 1 : 0;
    } else if (strcmp(key, "extra") == 0) {
        n->extra_loss_db = parse_double(val);
    }
    // nepoznati "key" su ignorirani
}

// parsira 1 liniju topologije
static void parse_line(Node* n, char* line_no_indent) {
    // Tokenizira po razmaku: TYPE key=val key=val ...
    char* save = NULL;
    char* tok = strtok_r(line_no_indent, " \t", &save);
    if (!tok) {
        die("Prazna linija topologije");
    }
    n->type = parse_type(tok);

    while ((tok = strtok_r(NULL, " \t", &save)) != NULL) {
        char* eq = strchr(tok, '=');
        if (!eq) continue; // ignorira neispravne tokene
        *eq = '\0';

        const char* key = tok;
        const char* val = eq + 1;
        apply_kv(n, key, val);
    }
}

// računa gubitak splittera
static double splitter_loss_db(int ratio) {
    if (ratio <= 1) return 0.0;
    // idealni split gubitak(10log10) + insertion loss - npr. 1:8  -> ~9 dB + insertion
    return 10.0 * log10((double)ratio) + SPLITTER_INS_DB;
}

// računa gubitke fizičkog optičkog linka 
static double node_link_loss_db(const Node* n) {
    double loss = 0.0;
    loss += n->len_km * ATTEN_DB_PER_KM;
    loss += (double)n->connectors * CONN_LOSS_DB;
    loss += (double)n->splices * SPLICE_LOSS_DB;

    if (n->type == NODE_SPLITTER) {
        loss += splitter_loss_db(n->splitter_ratio);
    }
    if (n->faulty) {
        loss += n->extra_loss_db;
    }
    return loss;
}

// inicijalizacija strukutre statistike
static SubtreeStats stats_init(void) {
    SubtreeStats s;

    s.ont_count = 0;
    s.ok_count = 0;
    s.fail_count = 0;
    s.down_count = 0;
    s.sum_rx = 0.0;
    s.sum_loss = 0.0;
    s.best_rx = -1e9;
    s.worst_rx =  1e9;
    return s;
}

// spaja statistiku djece u roditelja
static void stats_merge(SubtreeStats* a, const SubtreeStats* b) {
    a->ont_count += b->ont_count;
    a->ok_count += b->ok_count;
    a->fail_count += b->fail_count;
    a->down_count += b->down_count;
    a->sum_rx += b->sum_rx;
    a->sum_loss += b->sum_loss;
    if (b->best_rx > a->best_rx) {
        a->best_rx = b->best_rx;
    }
    if (b->worst_rx < a->worst_rx) {
        a->worst_rx = b->worst_rx;
    }
}

// dodaje putanje (OLT-SPLITTER-ONT)
static void path_append(char* path, size_t cap, const char* part) {
    if (path[0] != '\0') {
        strncat(path, "/", cap - strlen(path) - 1);
    }
    strncat(path, part, cap - strlen(path) - 1);
}

// rekurzivno prolazi topologiju od OLT-a prema ONT-ovima
static SubtreeStats walk_and_compute(
    const Node* n,
    double parent_tx_dbm,
    double rxmin_dbm,
    double acc_loss_db,
    double acc_dist_km,
    int down_flag,
    FILE* ont_csv,
    SplitterList* splitters,
    char* path,
    size_t path_cap
) {
    if (!n) return stats_init();

    // nasljeđuje OLT tx/rxmin ako je čvor OLT
    double tx_dbm = parent_tx_dbm;
    double my_rxmin = rxmin_dbm;

    if (n->type == NODE_OLT) {
        tx_dbm = n->olt_tx_dbm;
        my_rxmin = n->gpon_rxmin_dbm;
    }

    // za čvorove koji nisu OLT, dodaje njihov link loss u akumulirani loss
    double new_loss = acc_loss_db;
    double new_dist = acc_dist_km;
    int new_down = down_flag;

    if (n->type != NODE_OLT) {
        new_loss += node_link_loss_db(n);
        new_dist += n->len_km;
        if (n->faulty) {
            new_down = 1; // ako je neki element na putu faulty, svi ONT-ovi ispod se smatraju DOWN(idalje računa rx, ali je down)
        }
    }

    // ažurira path
    char old_path[512];
    strncpy(old_path, path, sizeof(old_path) - 1);
    old_path[sizeof(old_path) - 1] = '\0';

    char part[128];
    part[0] = '\0';
    if (n->type == NODE_OLT) {
        snprintf(part, sizeof(part), "OLT");
    } else if (n->type == NODE_SPLITTER) {
        if (n->name[0]) {
            snprintf(part, sizeof(part), "%s(1:%d)", n->name, n->splitter_ratio);
        } else {
            snprintf(part, sizeof(part), "S(1:%d)", n->splitter_ratio);
        }
    } else if (n->type == NODE_ONT) {
        snprintf(part, sizeof(part), "ONT#%d", n->ont_id);
    }

    if (part[0]) {
        path_append(path, path_cap, part);
    }

    SubtreeStats here = stats_init();

    if (n->type == NODE_ONT) {
        // racuna RX na ONT-u
        double rx_dbm = tx_dbm - new_loss;
        double margin = rx_dbm - my_rxmin;

        if (ont_results_count < MAX_ONT) {
            ont_results[ont_results_count].ont_id = n->ont_id;
            ont_results[ont_results_count].rx_dbm = rx_dbm;
            ont_results[ont_results_count].margin_db = margin;
            strncpy(ont_results[ont_results_count].path, path, 511);
            ont_results_count++;
        }

        const char* status;
        int ok = 0;
        if (new_down) {
            status = "DOWN";
        } else if (rx_dbm >= my_rxmin) {
            status = "OK";
            ok = 1;
        } else {
            status = "FAIL";
        }

        fprintf(ont_csv, "%d,%.4f,%.4f,%.4f,%.4f,%s,\"%s\"\n",
            n->ont_id, new_dist, new_loss, rx_dbm, margin, status, path);

        here.ont_count = 1;
        here.sum_rx = rx_dbm;
        here.sum_loss = new_loss;
        here.best_rx = rx_dbm;
        here.worst_rx = rx_dbm;

        if (new_down) {
            here.down_count = 1;
        } else if (ok) {
            here.ok_count = 1;
        } else {
            here.fail_count = 1;
        }

        // vracamo path
        strncpy(path, old_path, path_cap - 1);
        path[path_cap - 1] = '\0';
        return here;
    }

    // rekurzija djece
    const Node* child = n->child;
    while (child) {
        SubtreeStats cs = walk_and_compute(child, tx_dbm, my_rxmin, new_loss, new_dist, new_down, ont_csv, splitters, path, path_cap);
        stats_merge(&here, &cs);
        child = child->sibling;
    }

    // ako je čvor splitter, spremi podatke te grane (ONT-ovi su ispod)
    if (n->type == NODE_SPLITTER) {
        SplitterRecord rec;
        memset(&rec, 0, sizeof(rec));
        if (n->name[0]) {
            strncpy(rec.name, n->name, sizeof(rec.name) - 1);
        } else {
            strncpy(rec.name, "(unnamed)", sizeof(rec.name) - 1);
        }

        rec.ratio = n->splitter_ratio;
        rec.ont_count = here.ont_count;
        rec.ok_count = here.ok_count;
        rec.fail_count = here.fail_count;
        rec.down_count = here.down_count;
        rec.avg_rx = (here.ont_count > 0) ? (here.sum_rx / here.ont_count) : 0.0;
        rec.avg_loss = (here.ont_count > 0) ? (here.sum_loss / here.ont_count) : 0.0;
        rec.worst_rx = (here.ont_count > 0) ? here.worst_rx : 0.0;

        splitter_list_push(splitters, &rec);
    }

    // vracamo path
    strncpy(path, old_path, path_cap - 1);
    path[path_cap - 1] = '\0';
    return here;
}

// generiranje csv datoteke za splittere
static void write_splitter_csv(const char* filename, const SplitterList* sl) {
    FILE* f = fopen(filename, "w");
    if (!f) {
        die("Nemoguce je otvoriti splitter_results.csv za pisanje.");
    }

    fprintf(f, "name,ratio,ont_count,ok_count,fail_count,down_count,avg_rx_dbm,avg_loss_db,worst_rx_dbm\n");
    for (size_t i = 0; i < sl->n; i++) {
        const SplitterRecord* r = &sl->arr[i];
        fprintf(f, "\"%s\",%d,%d,%d,%d,%d,%.4f,%.4f,%.4f\n",
            r->name, r->ratio, r->ont_count, r->ok_count, r->fail_count, r->down_count, r->avg_rx, r->avg_loss, r->worst_rx);
    }
    fclose(f);
}

// ispis na konzolu
static void print_summary(const SubtreeStats* all, double tx, double rxmin) {
    printf("\n=== SUMMARY ===\n");
    printf("OLT TX: %.2f dBm | GPON RXmin: %.2f dBm\n", tx, rxmin);
    printf("ONT total: %d\n", all->ont_count);
    printf("OK:   %d\n", all->ok_count);
    printf("FAIL: %d\n", all->fail_count);
    printf("DOWN: %d\n", all->down_count);

    if (all->ont_count > 0) {
        printf("AVG RX:   %.2f dBm\n", all->sum_rx / all->ont_count);
        printf("AVG LOSS: %.2f dB\n", all->sum_loss / all->ont_count);
        printf("BEST RX:  %.2f dBm\n", all->best_rx);
        printf("WORST RX: %.2f dBm\n", all->worst_rx);
        printf("Note: Margin = RX - RXmin; PASS if RX >= RXmin and not DOWN.\n");
    }
}

// sortiranje ONT-ova uzlazno po optičkoj margini
int cmp_margin(const void* a, const void* b) {
    OntResult* x = (OntResult*)a;
    OntResult* y = (OntResult*)b;

    if (x->margin_db < y->margin_db) return -1;
    if (x->margin_db > y->margin_db) return 1;
    return 0;
}

// čitanje datoteke topologije
static Node* read_topology(const char* filename) {
    FILE* f = fopen(filename, "r");
    if (!f) {
        die("Nemoguce je otvoriti datoteku ftth_topology.txt");
    }
    Node* root = NULL;
    Node* stack[64] = {0};  //stack[depth] = zadnji cvor na depth

    char line[1024];
    int auto_ont_id = 1;

    while (fgets(line, sizeof(line), f)) {
        rtrim(line);
        char* p = ltrim(line);

        if (*p == '\0') continue;
        if (*p == '#') continue;

        int spaces = leading_spaces(line);
        if (spaces % 2 != 0) {
            die("Uvlaka mora biti paran broj razmaka");
        }

        int depth = spaces / 2;
        char tmp[1024];
        strncpy(tmp, p, sizeof(tmp) - 1);
        tmp[sizeof(tmp) - 1] = '\0';

        Node* n = node_new(NODE_ONT); // vrstu postavlja parse_line
        parse_line(n, tmp);

        if (n->type == NODE_ONT) {
            if (n->ont_id < 0) {
                n->ont_id = auto_ont_id++;
            }
        }

        if (depth == 0) {
            // OLT mora biti na samom vrhu
            if (n->type != NODE_OLT) {
                die("Najgornji cvor mora biti OLT!");
            }
            root = n;
            stack[0] = n;
        } else {
            Node* parent = stack[depth - 1];
            if (!parent) {
                die("Kriva identacija / Fali roditelj");
            }
            node_add_child(parent, n);
            stack[depth] = n;
        }

        // omogućava oslobadaje dubljeg satcka
        for (int i = depth + 1; i < 64; i++) {
            stack[i] = NULL;
        }
    }

    fclose(f);

    if (!root) {
        die("Nema OLT cvora u ftth_topology.txt");
    }

    return root;
}

void generate_report(const SubtreeStats* stats, Node* root) {
    FILE* f = fopen("report.txt", "w");
    if (!f) return;

    fprintf(f, "FTTH/GPON SIMULATION REPORT\n\n");
    fprintf(f, "OLT TX power: %.2f dBm\n", root->olt_tx_dbm);
    fprintf(f, "GPON RX minimum: %.2f dBm\n\n", root->gpon_rxmin_dbm);

    fprintf(f, "Total ONT count: %d\n", stats->ont_count);
    fprintf(f, "OK connections: %d\n", stats->ok_count);
    fprintf(f, "FAIL connections: %d\n", stats->fail_count);
    fprintf(f, "DOWN connections: %d\n\n", stats->down_count);

    fprintf(f, "Best RX power: %.2f dBm\n", stats->best_rx);
    fprintf(f, "Worst RX power: %.2f dBm\n", stats->worst_rx);
    fprintf(f, "Average RX power: %.2f dBm\n\n",
            stats->sum_rx / stats->ont_count);

    fprintf(f, "TOP %d worst ONT connections (by margin):\n", TOP_N);
    for (int i = 0; i < TOP_N && i < ont_results_count; i++) {
        fprintf(
            f,
            "ONT %d | margin = %.2f dB | RX = %.2f dBm | %s\n",
            ont_results[i].ont_id,
            ont_results[i].margin_db,
            ont_results[i].rx_dbm,
            ont_results[i].path
        );
    }

    fprintf(f, "\nZakljucak:\n");
    fprintf(f, "Veliki omjer dijeljenja(split ratio) i dugacke opticke udaljenosti smanjuju opticku marginu.\n");
    fprintf(f, "Kriticni dijelovi mreze bi se trebali pratiti i sanirti po potrebi!\n");

    fclose(f);
}

static void free_tree(Node* n) {
    if (!n) return;

    free_tree(n->child);
    free_tree(n->sibling);
    free(n);
}