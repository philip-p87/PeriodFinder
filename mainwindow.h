#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QSet>
#include <QMultiMap>

#define BUFFER_SIZE 128

#define FIND_PERIOD_MODIFIED_MODE 1
#define FIND_MAX_N_MODIFIED_MODE 2
#define FIND_SHORTEST_DEADLINE_MODIFIED_MODE 3      //modified mode = only n packets
#define FIND_PERIOD_OLD_MODE 4                      //old mode== many packets (2n+1)
#define FIND_MAX_N_OLD_MODE 5
#define FIND_SHORTEST_DEADLINE_OLD_MODE 6


struct result_t                     //for statistics
{
    unsigned int n;
    unsigned int period_buf[BUFFER_SIZE];          //holds values for periods (biggest period, lower index)
    unsigned int packet_size_list[BUFFER_SIZE];     //packet sizes from alternative config list
    unsigned int deadline_list[BUFFER_SIZE];        //deadlines from alternative config list
    unsigned int packet_number_list[BUFFER_SIZE];   //number of packets that have to be sent
    double deadline_list_normalized[BUFFER_SIZE];   //normalized deadlines from alternative config list
    double packet_duration_list[BUFFER_SIZE];       //packet durations from alternative config
    double packet_duration_list_normalized[BUFFER_SIZE];
    int sort_modus;                         //0 = none; 1 = shortest first; 2 = longest first
    double deadline_list_sorted[BUFFER_SIZE];
};

struct plot_t
{
    unsigned int step_size;
    unsigned int steps;
    unsigned int last_mode; //1 = DEEP Heur, 2 = UNUSED, 3 = DEEP analytic, 4 = EMAC
    unsigned int number_of_nodes[BUFFER_SIZE];
    unsigned int deadline[BUFFER_SIZE];
    unsigned int last_best_deadline[BUFFER_SIZE];
    unsigned long elapsed_time[BUFFER_SIZE];
};

struct EMAC_sequence_t
{
    unsigned int numbers[BUFFER_SIZE];
    //unsigned int sum_combinations[BUFFER_SIZE];
    unsigned int deadline;
    //int id;

    unsigned int number_of_sums;
    QList<unsigned int> sum_combinations_list;

    unsigned int n_conflict_with;
    QList<unsigned int> conflict_with_list;
};



namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

public slots:
    void saveText();
    void clearText();
    void updateSavePath();
    void aboutButtonPressed();

private slots:

    //DEEP heuristic
    bool calculatePeriodsDEEPH(unsigned int n);                //returns false if no periods are found
    int calculateMaxNPerD(int seed);                      //returns max number of nodes that fit in a specified deadline
    int calculateMaxNPerDDeepSearch(int seed, double d_max, double d_min, double* best_deadline);    //searches for highest N between d_min and d_max

    //DEEP analytic
    void calculatePeriodsDEEPA(unsigned int n);
    int calculateMaxNPerDDEEPA(int seed, unsigned int max_deadline, unsigned int packet_size_);

    //EMAC
    void EMAC_findPeriodButtonPressed();

    //Plot functions
    void plotMaxNPerDeadlineDEEPHPressed();
    void plotMaxNPerDeadlineDEEPAPressed();
    void plotMaxNPerDeadlineEMACPressed();
    void plotSaveFromBufToFileButtonPressed();
    //new 28 Dec 21
    void plotMaxNPerDeadlinePKPressed();
    int pkCalculateMaxN(int seed, unsigned int max_deadline, unsigned int packet_size_);
    int pkMaxDeadline(int n, unsigned int packet_size_);

    //create period list functions (work in progress)
    void createPeriodListDEEPHeuristicButtonPressed();
    void createPeriodListDEEPAnalyticButtonPressed();
    void createPeriodListEMACButtonPressed();


private:
    Ui::MainWindow *ui;

    inline bool check_if_period_is_valid(unsigned int period_buf_index, unsigned int n);
    void calculatePacketNumbers(unsigned int n);
    void updatePacketDeadlineList();                      //read ui input forms


    //EMAC
    bool EMAC_calculatePeriods(unsigned int deadline_given, unsigned int n_given, bool debug_output);          //returns max number of nodes. All other values can be accessed through global variables
    void combinationUtil(unsigned int arr[], int n, int r, int index, int data[], int i);       //Program to print all combination of size r in an array of size n
    bool areDisjoint(unsigned int index_1, unsigned int index_2);                               //checks if two arrays have values in common
    QSet<unsigned int> sum_combinations;                            //string period combinations. temporarily used for gnerating sums
    QMultiMap<unsigned int, unsigned int> n_conflict_with_map;      //QMap <n_conflict_with, sequence_id>  stores the number of conflicts and id
    QMultiMap<unsigned int, unsigned int> deadline_map;             //QMap <deadline, sequence_id>
    EMAC_sequence_t a_sequence;
    QList<EMAC_sequence_t> sequence_list;
    //EMAC_sequence_t selected_sequence_list[100];
    QList<EMAC_sequence_t> selected_sequence_list2;
    QList<EMAC_sequence_t> shortest_selected_sequence_list;


    //helping functions
    int lcm(int a, int b);
    int gcd(int a, int b);

    //variables
    unsigned int period_buf[BUFFER_SIZE];          //holds values for periods (biggest period, lower index)
    unsigned int period_buf_temp_[BUFFER_SIZE];    //temporary values
    //note that the following values are stored in an array rather than a single value
    //this is for testing (in particular different deadlines, etc.)
    unsigned int packet_size_list[BUFFER_SIZE];     //packet sizes from alternative config list
    unsigned int deadline_list[BUFFER_SIZE];        //deadlines from alternative config list
    unsigned int packet_number_list[BUFFER_SIZE];   //number of packets that have to be sent
    double deadline_list_normalized[BUFFER_SIZE];   //normalized deadlines from alternative config list
    double packet_duration_list[BUFFER_SIZE];       //packet durations from alternative config
    double packet_duration_list_normalized[BUFFER_SIZE];

    //storages for plotting functions
    plot_t plot_storage;

    //storage container for saving results to text
    unsigned int storage_index;
    result_t storage[BUFFER_SIZE];

    //misc
    void setBusyState(bool state); //set indicator in gui that shows if program is busy or not

};

#endif // MAINWINDOW_H
