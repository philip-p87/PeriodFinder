#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <stdlib.h>
#include <QtGui>
#include <QFileDialog>
#include <QTextCursor>
#include <QElapsedTimer>
#include <QMessageBox>

#include <math.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <algorithm>


QString save_location = "C:/Users/Philip/Desktop/output.csv";


MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    setWindowTitle("Reliable Unidirectional Communication -- Period Finder");
    ui->pathEdit->setText(save_location);

    connect(ui->clearButton, SIGNAL(clicked()), this, SLOT(clearText()));
    connect(ui->saveButton, SIGNAL(clicked()), this, SLOT(saveText()));
    connect(ui->pathEdit, SIGNAL(editingFinished()), this, SLOT(updateSavePath()));
    connect(ui->aboutButton, SIGNAL(clicked()), this, SLOT(aboutButtonPressed()));

    connect(ui->plotMaxNPerDeadlineDEEPH, SIGNAL(clicked()), this, SLOT(plotMaxNPerDeadlineDEEPHPressed()));
    connect(ui->plotMaxNPerDeadlineDEEPA, SIGNAL(clicked()), this, SLOT(plotMaxNPerDeadlineDEEPAPressed()));
    connect(ui->plotMaxNPerDeadlineEMAC, SIGNAL(clicked()), this, SLOT(plotMaxNPerDeadlineEMACPressed()));
    connect(ui->plotMaxNPerDeadlinePK, SIGNAL(clicked()), this, SLOT(plotMaxNPerDeadlinePKPressed()));

    connect(ui->createPeriodListDEEPHeuristicButton, SIGNAL(clicked()), this, SLOT(createPeriodListDEEPHeuristicButtonPressed()));
    connect(ui->createPeriodListDEEPAnalyticButton, SIGNAL(clicked()), this, SLOT(createPeriodListDEEPAnalyticButtonPressed()));
    connect(ui->createPeriodListEMACButton, SIGNAL(clicked()), this, SLOT(createPeriodListEMACButtonPressed()));


    setBusyState(false);
}

MainWindow::~MainWindow()
{
    //free(buffer);
    delete ui;
}



//this function searches for periods (aka inter-packet times). Returns true, if successful. Results are stored in period_buf[].
//n = number of nodes
bool MainWindow::calculatePeriodsDEEPH(unsigned int n)
{
    bool found_period_flag;
    bool period_check;  //for debugging
    period_buf[0] = (deadline_list_normalized[0]-packet_size_list[0]) / (packet_number_list[0]-1);
    period_buf_temp_[0] = period_buf[0];

    for(unsigned int x=1; x<n; x++)    //TODO: corrected x<0 to x<=0 -> check results on other functions <- THIS WAS WRONG, changed it back and it works now
    {
        found_period_flag = false;
        for (int buf = (deadline_list_normalized[x]-packet_size_list[x])/(packet_number_list[x]-1); buf > 1; buf--)
        {
            period_buf_temp_[x] = buf;
            period_check = check_if_period_is_valid(x, n);
            if(period_check == true)
            {
                found_period_flag = true;
                period_buf[x] = buf;
                break;
            }
        }
        if(found_period_flag == false)
        {
            //did not find all period -> fill period storage with zeros
            for(unsigned int i=x; i<=n; i++)
                period_buf[i] = 0;
            //break;
            return false;
        }
    }
    return true;
}

inline bool MainWindow::check_if_period_is_valid(unsigned int period_buf_index, unsigned int n)
{
    for(unsigned int j=0;j<period_buf_index; j++)
    {
        int p_j = period_buf_temp_[j];
        int p_i = period_buf_temp_[period_buf_index];
        int l_j = packet_size_list[j];
        int l_i = packet_size_list[period_buf_index];

        if(p_j < p_i) //vorherige Periode war kleiner
        {
            for(unsigned int k = 1; k < n; k++)
            {
                double temp = (k * p_i) % p_j;
                if((temp < (l_i + l_j)) || (p_j - temp < (l_i + l_j)))
                {
                    return false;
                }
                //else
                //    continue;
            }
        }
        else //vorherige Periode war größer
        {
            for(unsigned int k = 1; k < n; k++)
            {
                double temp = (k * p_j) % p_i;
                if((temp < (l_i + l_j)) || (p_i - temp < (l_i + l_j)))
                {
                    return false;
                }
            }
        }


    }
    return true;
}




//read deadline, packet length etc. from gui
void MainWindow::updatePacketDeadlineList()
{
    for(int i=0; i<BUFFER_SIZE; i++)
    {
        packet_size_list[i] = ui->packetSizeEdit->text().toUInt();
        deadline_list[i] = ui->deadlineEdit->text().toUInt();
        packet_duration_list[i] = (double) packet_size_list[i] *8 / ui->speedEdit->text().toDouble();
        deadline_list_normalized[i] = deadline_list[i] /(8 / ui->speedEdit->text().toDouble());
    }
}

//returns number of nodes that fit into specified deadline (which must be separately be filled in deadline_list_normalized,etc)
//seed gives the starting value for search --> used for accelerating the search
int MainWindow::calculateMaxNPerD(int seed =0)
{
    bool find_period_flag;
    unsigned int n;

    //finding periods
    if(seed == 0)
        n = 3;
    else
        n = seed;


    for(; n < BUFFER_SIZE; n++)
    {
        //calculate packet numbers
        calculatePacketNumbers(n);

        //calculate periods
        find_period_flag = calculatePeriodsDEEPH(n);
        if(find_period_flag == false)
        {
            //calculate last valid n
            n = n - 1;
            calculatePacketNumbers(n);

            find_period_flag = calculatePeriodsDEEPH(n);

            if(find_period_flag == false)
                return 0;
            return n;
        }


        //update GUI to prevent GUI freezing
        QCoreApplication::processEvents(QEventLoop::AllEvents, 10);
    }
    return 0;//found n greater than BUFF size
}

int MainWindow::calculateMaxNPerDDEEPA(int seed = 3, unsigned int max_deadline = 0, unsigned int packet_size_ = 0)
{
    unsigned int n;
    unsigned int deadline;
    unsigned int packet_size;
    if(packet_size_ == 0)
        packet_size = ui->packetSizeEdit->text().toInt();
    else
        packet_size = packet_size_;
    double transmission_duration = 8 / ui->speedEdit->text().toDouble();    //duration of time base (8bit)

    //convert to normalized
    max_deadline = max_deadline / transmission_duration;

    for(n = seed; n < BUFFER_SIZE; n++)
    {
        //          -period------------  -longest period- -l_max--------- -period count-
        deadline = ((2*(n-2)*(n-1) + 2) + (n-1)*2)        * packet_size * (n-1);

        if(deadline > max_deadline)
        {
            //calculate last valid n
            n = n - 1;
            return n;
        }


        //update GUI to prevent GUI freezing
        QCoreApplication::processEvents(QEventLoop::AllEvents, 10);
    }
    return 0;//found n greater than BUFF size
}


int MainWindow::calculateMaxNPerDDeepSearch(int seed = 3, double d_max = 100, double d_min = 0, double *best_deadline = NULL)
{
    bool find_period_flag;
    int n;
    int n_highest = seed;
    int d_min_normalized, d_max_normalized;
    double time_base = 8 / ui->speedEdit->text().toDouble();    //duration of time base (8bit)

    d_min_normalized = d_min / time_base;
    d_max_normalized = d_max / time_base;

    if(d_min < 0)
        d_min = 0;


    //for(double deadline = d_min_normalized; deadline < d_max_normalized; deadline += 1)
    for(double deadline = d_min_normalized; deadline < d_max_normalized; deadline += ui->packetSizeEdit->text().toDouble()) //test debug
    {
        //update deadlines
        for(int i=0; i< BUFFER_SIZE; i++)
        {
            deadline_list_normalized[i] = deadline;
        }


        //search for max n
        for(n = seed; n < BUFFER_SIZE; n++)
        {
            //calculate packet numbers for normal and modified journal alg
            calculatePacketNumbers(n);

            //calculate periods
            find_period_flag = calculatePeriodsDEEPH(n);
            if(find_period_flag == false)
            {
                n--;
                if(n >= n_highest)
                {
                    n_highest = n;
                    //if(&best_deadline != NULL)
                    *best_deadline = (double)deadline;
                }
                break;
            }

            //update GUI to prevent GUI freezing
            QCoreApplication::processEvents(QEventLoop::AllEvents, 10);
        }
    }

*best_deadline = (*best_deadline) * time_base; //conversion to ms
return n_highest;
}


//calculates packet numbers for DEEP
void MainWindow::calculatePacketNumbers(unsigned int n)
{
    unsigned int packet_number_buf;
    for(unsigned int i=0; i<n; i++)
    {
        packet_number_buf = 0;
        for(unsigned int k=0; k<n; k++)
        {
            if(k==i)
                continue;

            packet_number_buf += deadline_list[i]/deadline_list[k]; //gamma formel für journal old (ohne inter-sequence pause)
        }
        packet_number_list[i] = n;
    }
}

void MainWindow::calculatePeriodsDEEPA(unsigned int n)
{
    unsigned int l_max;

    packet_number_list[0] = n;
    packet_size_list[0] = ui->packetSizeEdit->text().toInt();
    l_max = packet_size_list[0];
    period_buf[0] = (2*(n-2)*(n-1) + 2) * l_max;
    for(unsigned int i=1; i<n; i++)
    {
        packet_size_list[i] = packet_size_list[0];
        packet_number_list[i] = packet_number_list[0];
        deadline_list[i] = deadline_list[0];

        period_buf[i] = period_buf[i-1] + 2*l_max;
    }
}





////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////// EMAC
///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void MainWindow::EMAC_findPeriodButtonPressed()
{
    EMAC_calculatePeriods(0,0,false);
    ui->logEdit->appendPlainText(QString("\n Number of sequence_list: %1      Number of selected_sequence_list2:%2\n")
                                 .arg(sequence_list.count()).arg(selected_sequence_list2.count()));
}



//returns max numer of nodes that can fit in deadline
//input: deadline= if 0 then uses gui deadline
//input: debug_output=true print debug messages to gui
bool MainWindow::EMAC_calculatePeriods(unsigned int deadline_given = 0, unsigned int n_given = 0, bool debug_output = false)
{
    int n, counter = 0;
    unsigned int buffer;
    double deadline, time_base;
    bool success = false;       //this value is returned: true = enough periods were found

    time_base = 1 / ui->speedEdit->text().toDouble();       //transmission duration of 1byte //TODO change back

    if(deadline_given > 0)
        deadline = deadline_given;
    else
        deadline = ui->deadlineEdit->text().toInt() / (time_base * ui->packetSizeEdit->text().toDouble());

    if(n_given != 0)
        n = n_given;
    else
        n = 2;//ui->minBox->text().toInt();


    //clear old data////////////////////////////////////////////////////////////////////////////////
    shortest_selected_sequence_list.clear();
    deadline_map.clear();
    n_conflict_with_map.clear();
    for (int i = 0; i < sequence_list.size(); ++i)
    {
        //sequence_list[i].numbers.clear();
        sequence_list[i].sum_combinations_list.clear();
        sequence_list[i].conflict_with_list.clear();

    }
    for (int i = 0; i < selected_sequence_list2.size(); ++i)
    {
        //selected_sequence_list2[i].numbers.clear();
        selected_sequence_list2[i].sum_combinations_list.clear();
        selected_sequence_list2[i].conflict_with_list.clear();

    }
    //a_sequence.numbers.clear();
    a_sequence.sum_combinations_list.clear();
    a_sequence.conflict_with_list.clear();
    sequence_list.clear();
    selected_sequence_list2.clear();


////////////////////////////////////////////////////////////////////////
////////Part 1: Create Sequences
///////////////////////////////////
    //little performance tweak
    int factor_maximum = 150;
    if(n>14)
        factor_maximum += (n-14)*10;

    //EMAC original
    int NREPL = n-1;
    //int nsequences = 0;
    for (int i=2;i<=factor_maximum;i+=2)
    {
      //sequences[nsequences].multiplier = i;
      for (int k=0;k<NREPL;k++)
      {
        //sequences[nsequences].thenumbers[k] = sequences[nsequences].multiplier;
          a_sequence.numbers[k] = i;
      }

      //if (getsum_in_create_seq(nsequences)+1<=MAXZ)
        //nsequences = nsequences + 1;
      //calculate sum (i.e. sequence length)
      buffer = 0;
      for(int j = 0; j<n-1; j++)
      {
          buffer += a_sequence.numbers[j];
      }
      buffer++; //last packet
      a_sequence.deadline = buffer;

      //sequence valid? If yes, copy sequence, else skip iteration
      if(buffer <= deadline)
      {
          sequence_list.append(a_sequence);
          counter++;
      }
    }


////////////////////////////////////////////////////////////////////////
////////Part 2: Calculate all possible period combinations
///////////////////////////////////



    int data[10000];

    for(int index = 0; index < counter; index++)
    {
        //different combination pairs, for example k=2 will give all combinations of 2 elements
        for(int k = 1; k <= n; k++)
        {
            // Print all combination using temprary array 'data[]', Will print results in QList
            combinationUtil(sequence_list[index].numbers, n-1, k, 0, data, 0);
        }

        //convert to QList for sorting
        sequence_list[index].number_of_sums = sum_combinations.count();
        sequence_list[index].sum_combinations_list = sum_combinations.toList();
        qSort(sequence_list[index].sum_combinations_list);
        sum_combinations.clear();
    }





////////////////////////////////////////////////////////////////////////
////////Part 3: Count number of collisions
///////////////////////////////////



    for(int index = 0; index < counter; index++)
    {
        //initialize with 0
        sequence_list[index].n_conflict_with = 0;

        //count all combinations with all other sequences
        for(int j = 0; j < counter; j++)
        {
            if(j == index)
                continue;

            if(areDisjoint(index, j) == false)  //false if both arrays have at least an element in common (collision)
            {
                sequence_list[index].n_conflict_with++;
                sequence_list[index].conflict_with_list.append(j);
            }
        }

        //store n_conflict_with and id for PART 4
        n_conflict_with_map.insert(sequence_list.at(index).n_conflict_with, (unsigned int) index);
    }

////////////////////////////////////////////////////////////////////////
////////Part 4: Select suitable sequences and fill them in select_sequence_list
///////////////////////////////////



    bool sequence_ok_flag;      //false: collisions, true: no collisions
    int selected_counter = 0;

    //The Multimap contains Key: n_conflicts_with and value: sequence_list id.
    //   Key was just used for sorting (now obsolete), value will be used for accessing the sequence
    QMultiMap<unsigned int, unsigned int>::const_iterator selected_sequence_iterator = n_conflict_with_map.constBegin();
    while (selected_sequence_iterator != n_conflict_with_map.constEnd())
    {
        sequence_ok_flag = true;

        //iterate through all selected_sequences and check if there is one that collides with the current sequence from sequence_list
        for (int j = 0; j < selected_sequence_list2.size(); ++j)
        {
            if(selected_sequence_list2.at(j).conflict_with_list.contains(selected_sequence_iterator.value()))
            {
                sequence_ok_flag = false;       //we found a selected_sequence that collides with our current iteration -> set flag to false
                break;
            }
        }

        if(sequence_ok_flag)
        {
            //sequence was ok -> copy it to selected_sequence_list
            selected_sequence_list2.append(sequence_list[selected_sequence_iterator.value()]);

            deadline_map.insert(sequence_list[selected_sequence_iterator.value()].deadline, selected_counter);
            selected_counter++;
        }

        ++selected_sequence_iterator;
    }


    /*if(debug_output)
    {
        //Debug
        ui->logEdit->insertPlainText(QString("EMAC Paper debug output\n"));
        ui->logEdit->insertPlainText(QString("n= %1   deadline = %2 \n").arg(n).arg(deadline));
        for (int i = 0; i < selected_sequence_list2.size(); ++i)
        {
            ui->logEdit->insertPlainText(QString("selected_sequence_list2[%1].numbers = ").arg(i));
            for(int j = 0;  j < n-1 ; j++)
            {
                ui->logEdit->insertPlainText(QString("%1, ").arg(selected_sequence_list2.at(i).numbers[j]));
            }
            ui->logEdit->insertPlainText(QString("    deadline %1, ").arg(selected_sequence_list2.at(i).deadline));
            ui->logEdit->insertPlainText(QString("\n"));

        }

        ui->logEdit->insertPlainText(QString("-----------------------------------------------\n"));
    }*/

////////////////////////////////////////////////////////////////////////
////////Part 5: Select best sequences from select_sequence_list
///////////////////////////////////



    int shortest_sequence_counter = 0;
    QMultiMap<unsigned int, unsigned int>::const_iterator deadline_iterator = deadline_map.constBegin();
    while (deadline_iterator != deadline_map.constEnd())
    {
        shortest_sequence_counter++;
        shortest_selected_sequence_list.append(selected_sequence_list2.at(deadline_iterator.value()));
        ++deadline_iterator;

        if(shortest_sequence_counter >= n)
        {
            success = true;
            break;
        }
    }


////////////////////////////////////////////////////////////////////////
////////Part 6: Debug output
///////////////////////////////////

    if(debug_output)
    {
        //Debug
        ui->logEdit->insertPlainText(QString("Selected sequences: deadline %1, n selected %2\n")
                                     .arg(deadline).arg(n));
        //ui->logEdit->insertPlainText(QString("n= %1   deadline = %2 \n").arg(n).arg(deadline));
        for (int i = 0; i < shortest_selected_sequence_list.size(); ++i)
        {
            ui->logEdit->insertPlainText(QString("shortest_selected_sequence_list[%1].numbers = ").arg(i));
            for(int j = 0;  j < n-1 ; j++)
            {
                ui->logEdit->insertPlainText(QString("%1, ").arg(shortest_selected_sequence_list.at(i).numbers[j]));
            }
            ui->logEdit->insertPlainText(QString("    deadline %1").arg(shortest_selected_sequence_list.at(i).deadline));
            ui->logEdit->insertPlainText(QString("\n"));


            //copy to buffer
            //period_buf[n-2] = shortest_selected_sequence_list.at(i).numbers[0];
            //ui->logEdit->insertPlainText(QString("period_buf[n-2] =%1 \n").arg(period_buf[n-2]));

        }

        ui->logEdit->insertPlainText(QString("-----------------------------------------------\n"));
    }

    /*//debug output
    ui->logEdit->insertPlainText(QString("EMAC Paper debug output\n"));
    ui->logEdit->insertPlainText(QString("n= %1   deadline = %2 \n").arg(n).arg(deadline));

    for(int i=0; i<counter; i++)
    {
        //ui->logEdit->insertPlainText(QString("sequence_list[%1].numbers    \n").arg(counter));

        //periods
        ui->logEdit->insertPlainText(QString("sequence_list[%1].numbers = ").arg(i));
        for(int j=0;j<n-1;j++)
        {
            ui->logEdit->insertPlainText(QString("%1, ").arg(sequence_list[i].numbers[j]));
        }
        ui->logEdit->insertPlainText(QString("    deadline %1, ").arg(sequence_list[i].deadline));
        ui->logEdit->insertPlainText(QString("\n"));


        //sums
        ui->logEdit->insertPlainText(QString("sequence_list[%1].sum_combinations_list = ").arg(i));
        for(unsigned int j=0;j<sequence_list[i].number_of_sums;j++)
        {
            ui->logEdit->insertPlainText(QString("%1, ").arg(sequence_list[i].sum_combinations_list.at(j)));
        }
        ui->logEdit->insertPlainText(QString("\n"));

        //number of conflicts
        ui->logEdit->insertPlainText(QString("sequence_list[%1].n_conflict_with = %2\n").arg(i).arg(sequence_list[i].n_conflict_with));

        //list of conflicts
        ui->logEdit->insertPlainText(QString("sequence_list[%1].conflict_with_list = ").arg(i));
        for(unsigned int j=0;j<sequence_list[i].n_conflict_with;j++)
        {
            ui->logEdit->insertPlainText(QString("%1, ").arg(sequence_list[i].conflict_with_list[j]));
        }
        ui->logEdit->insertPlainText(QString("\n\n"));

    }

    ui->logEdit->insertPlainText(QString("-----------------------------------------------\n"));*/
    return success;
}



/* arr[]  ---> Input Array
   n      ---> Size of input array
   r      ---> Size of a combination to be printed
   index  ---> Current index in data[]
   data[] ---> Temporary array to store current combination
   i      ---> index of current element in arr[]     */
void MainWindow::combinationUtil(unsigned int arr[], int n, int r, int index, int data[], int i)
{
    static unsigned int buf;
    // Current cobination is ready, print it
    if (index == r)
    {
        buf = 0;
        for (int j=0; j<r; j++)
        {
            buf += data[j];
        }
        sum_combinations.insert(buf);
        //printf("%d ",buf);
        //printf("\n");
        return;
    }

    // When no more elements are there to put in data[]
    if (i >= n)
        return;

    // current is included, put next at next location
    data[index] = arr[i];
    combinationUtil(arr, n, r, index+1, data, i+1);

    // current is excluded, replace it with next (Note that
    // i+1 is passed, but index is not changed)
    combinationUtil(arr, n, r, index, data, i+1);
}


// Returns true if set1[] and set2[] are disjoint, else false
bool MainWindow::areDisjoint(unsigned int index_1, unsigned int index_2)
{
    for (int i=0; i<sequence_list[index_2].sum_combinations_list.count(); i++)
        if (sequence_list[index_1].sum_combinations_list.contains(sequence_list[index_2].sum_combinations_list[i]))
            return false;

    return true;
}


int MainWindow::lcm(int a, int b)
{
  return (a*b)/gcd(a,b);
}

int MainWindow::gcd(int a, int b)
{
  if (b == 0)
    return a;
  else
    return gcd(b, a%b);
}





///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Plotfunktionen
///
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



void MainWindow::plotMaxNPerDeadlineDEEPHPressed()
{
    double min_deadline, max_deadline, step_size, deadline;
    int n = 0, n_old = 3, n_highest = 0;
    unsigned int k=0;
    long time_elapsed;
    QElapsedTimer timer;

    double comp_time[BUFFER_SIZE];
    int comp_time_n[BUFFER_SIZE];   //stores n for comp_time[]
    int comp_time_index = 0;

    setBusyState(true);

    //get information
    updatePacketDeadlineList();
    double time_base = 8 / ui->speedEdit->text().toDouble();
    //std::cout << "time_base: " << time_base << std::endl;

    //output csv ---------------------------------------------------------------------------------------
    std::ofstream myfile;
    QString message;
    myfile.open (save_location.toUtf8().constData(), std::ios::trunc);
    if(myfile.is_open() == false)
    {
        ui->logEdit->insertPlainText(QString("cannot open output log file!\n"));
        setBusyState(false);
        return;
    }
    min_deadline = ui->plotMaxNMinDEdit->text().toDouble();
    max_deadline = ui->plotMaxNMaxDEdit->text().toDouble();
    step_size = ui->plotMaxNStepSizeEdit->text().toDouble();

    plot_storage.step_size = step_size;

    if(step_size > min_deadline)
        deadline = 0;
    else
        deadline = min_deadline;

    if(ui->plotDeepSearchBox->isChecked()==false)
        message += QString("plotMaxNPerDeadlinePressed deep:%4, from d=%1 to d=%2, stepsize %3, timebase %5\n").arg(min_deadline)
            .arg(max_deadline).arg(step_size).arg(ui->plotDeepSearchBox->isChecked()).arg(time_base);
    else
        message += QString("plotMaxNPerDeadlinePressed deep:%4, from d=%1 to d=%2, stepsize %3, timebase %5, search delta:%6\%\n").arg(min_deadline)
            .arg(max_deadline).arg(step_size).arg(ui->plotDeepSearchBox->isChecked()).arg(time_base).arg(ui->searchDeltaEdit->text().toDouble());

    message += "\n";
    myfile << message.toStdString().data();
    ui->logEdit->appendPlainText(message);

    if(ui->plotDeepSearchBox->isChecked()==false)
        myfile << "deadline[ms]; Numbers of nodes; calculation time [ms]; periods [ms];\n";
    else
        myfile << "deadline[ms]; Numbers of nodes; calculation time [ms]; optimal deadline[ms]; delta[%]; periods [ms];\n";

    myfile.close();



    for(; deadline <= max_deadline; deadline = deadline + step_size, k++)
    {
        timer.restart();

        if(ui->plotDeepSearchBox->isChecked() == false)
        {
            //update input for calculateMaxNPerD()
            for(unsigned int i = 0; i < BUFFER_SIZE; i++)
            {
                deadline_list_normalized[i] = (int) (deadline / time_base);
            }

            //little optimization
            if(n > 10)
                n = calculateMaxNPerD(n_old-10);
            else
                n = calculateMaxNPerD(n);

            //TODO

            time_elapsed = timer.elapsed();

            if(n_old < n)
            {
                n_highest = n;
                n_old = n;
                comp_time[comp_time_index] = time_elapsed;
                comp_time_n[comp_time_index] = n_highest;
                comp_time_index++;
                //std::cout << "n_highest: " << n_highest << "   time_elapsed: " << time_elapsed << std::endl;
            }
            calculatePacketNumbers(n);
            calculatePeriodsDEEPH(n);
            ui->logEdit->insertPlainText(QString("deadline = %1 ms  max n = %2\n").arg(deadline).arg(n));
            plot_storage.number_of_nodes[k] = n;
            plot_storage.deadline[k] = (unsigned int) deadline;
            plot_storage.elapsed_time[k] = (int)timer.elapsed();

            myfile.open (save_location.toUtf8().constData(), std::fstream::out | std::fstream::app);
            if(myfile.is_open() == false)
            {
                ui->logEdit->insertPlainText(QString("cannot open output log file!\n"));
                setBusyState(false);
                return;
            }
            myfile << deadline << ";" << n << ";" << time_elapsed << ";" ;
            for(int i = 0; i<n; i++)
            {
                myfile << period_buf[i] * time_base << ";";
            }
            myfile << "\n";
            myfile.close();

        }

        //deep search
        else
        {
            if (deadline == max_deadline) break;
            double best_deadline = 0;       //deadline, at which the optimum n was found
            double search_delta = ui->searchDeltaEdit->text().toDouble()/100;
            double deadline_search_upper_bound;//, deadline_search_lower_bound;
            static double last_best_deadline;

            deadline_search_upper_bound = deadline + step_size;
            //n = calculateMaxNPerDDeepSearch(n_old, deadline + step_size, deadline, &best_deadline);
            n = calculateMaxNPerDDeepSearch(n_old, deadline_search_upper_bound, deadline_search_upper_bound *(1-search_delta), &best_deadline);

            time_elapsed = timer.elapsed();

            if(n<n_old) //new n smaller than previous value?
            {
                n = n_old;
                best_deadline = last_best_deadline;
            }

            last_best_deadline = best_deadline;
            ui->logEdit->insertPlainText(QString("deadline = %1 ms  max n = %2\n").arg(deadline + step_size).arg(n));
            plot_storage.number_of_nodes[k] = n;
            plot_storage.deadline[k] = (unsigned int) (deadline) + (unsigned int) (step_size);
            plot_storage.deadline[k] = last_best_deadline;
            plot_storage.elapsed_time[k] = (int)timer.elapsed();

            //save computation times
            if(n > n_highest)
            {
                n_highest = n;
                comp_time[comp_time_index] = time_elapsed;
                comp_time_n[comp_time_index] = n_highest;
                comp_time_index++;
            }

            //calculate periods (to save in file)
            for(unsigned int i = 0; i < BUFFER_SIZE; i++)
            {
                deadline_list_normalized[i] = best_deadline / time_base;
            }
            calculatePacketNumbers(n);
            calculatePeriodsDEEPH(n);

            //save file
            myfile.open (save_location.toUtf8().constData(), std::fstream::out | std::fstream::app);
            if(myfile.is_open() == false)
            {
                ui->logEdit->insertPlainText(QString("cannot open output log file!\n"));
                setBusyState(false);
                return;
            }

            double delta = (deadline + step_size - best_deadline)/(deadline + step_size) *100;
            myfile << deadline + step_size<< ";" << n << ";" << time_elapsed << ";" << best_deadline << ";" << delta << ";";
            for(int i = 0; i<n; i++)
            {
                myfile << period_buf[i] * time_base << ";";
            }
            myfile << "\n";
            myfile.close();
        }


    }
    plot_storage.steps = k;
    plot_storage.last_mode = 1;


    //save comp time in csv
    myfile.open (save_location.toUtf8().constData(), std::fstream::out | std::fstream::app);
    myfile << "\n\n" << "n; computation time;\n";
    for(int i=0; i<comp_time_index; i++)
    {
        myfile << comp_time_n[i] << ";" << comp_time[i] << ";\n";
    }
    myfile.close();

    ui->logEdit->appendPlainText(QString("finished!\n"));
    myfile.close();

    setBusyState(false);
}

void MainWindow::plotMaxNPerDeadlineDEEPAPressed()
{
    double min_deadline, max_deadline, step_size, deadline;
    int n, n_old = 2;
    unsigned int k=0;
    QElapsedTimer timer;
    long int time_elapsed;
    double comp_time[BUFFER_SIZE];
    int comp_time_n[BUFFER_SIZE];   //stores n for comp_time[]
    int comp_time_index = 0;
    QString message;

    setBusyState(true);

    std::ofstream myfile;
    myfile.open (save_location.toUtf8().constData(), std::ios::trunc);
    if(myfile.is_open() == false)
    {
        ui->logEdit->insertPlainText(QString("cannot open output log file!\n"));
        setBusyState(false);
        return;
    }
    min_deadline = ui->plotMaxNMinDEdit->text().toDouble();
    max_deadline = ui->plotMaxNMaxDEdit->text().toDouble();
    step_size = ui->plotMaxNStepSizeEdit->text().toDouble();

    if(step_size > min_deadline)
        deadline = 0;
    else
        deadline = min_deadline;


    message += QString("plotMaxNPerDeadlineDEEPAPressed from d=%1 to d=%2, stepsize %3\n").arg(min_deadline)
            .arg(max_deadline).arg(step_size);

    message += "\n";

    ui->logEdit->appendPlainText(message);
    myfile << message.toStdString().data();
    myfile << "deadline[ms]; Numbers of nodes; calculation time [ms]; periods\n";

    //get information
    updatePacketDeadlineList();
    double transmission_duration = 8 / ui->speedEdit->text().toDouble();    //duration of time base (8bit)

    for(; deadline <= max_deadline; deadline = deadline + step_size, k++)
    {
        timer.restart();

        //update input for calculateMaxNPerD()
        for(unsigned int i = 0; i < BUFFER_SIZE; i++)
        {
            deadline_list_normalized[i] = (int) (deadline / transmission_duration);
        }

        n = calculateMaxNPerDDEEPA(n_old, deadline);

        time_elapsed = timer.elapsed();
        if(n > n_old) //update compute time buf only if n increased
        {
            n_old = n;
            comp_time[comp_time_index] = time_elapsed;
            comp_time_n[comp_time_index] = n_old;
            comp_time_index++;
            //std::cout << "n_highest: " << n_highest << "   time_elapsed: " << time_elapsed << std::endl;
        }

        ui->logEdit->insertPlainText(QString("deadline = %1 ms  max n = %2\n").arg(deadline).arg(n));
        myfile << deadline << ";" << n << ";" << time_elapsed <<  ";";
        //append periods
        calculatePeriodsDEEPA(n);
        for(int i=0; i<n; i++) {
            myfile << period_buf[i] << ";";
        }
        myfile  <<  "\n";

    }

    //save comp time in csv
    myfile << "\n\n" << "n; computation time;\n";
    for(int i=0; i<comp_time_index; i++)
    {
        myfile << comp_time_n[i] << ";" << comp_time[i] << ";\n";
    }
    myfile.close();
    ui->logEdit->appendPlainText("finished!\n");

    setBusyState(false);
}

//lmax = packetSizeEdit
void MainWindow::plotMaxNPerDeadlineEMACPressed()
{
    double min_deadline, max_deadline, step_size, deadline;
    unsigned int k=0, n, nmax=0;
    double transmission_duration;
    QElapsedTimer timer;
    double comp_time[BUFFER_SIZE];
    int comp_time_n[BUFFER_SIZE];   //stores n for comp_time[]
    int comp_time_index = 0;
    long double time_elapsed;

    setBusyState(true);

    transmission_duration = ui->packetSizeEdit->text().toDouble() * 8 / ui->speedEdit->text().toDouble();
    min_deadline = ui->plotMaxNMinDEdit->text().toDouble();
    max_deadline = ui->plotMaxNMaxDEdit->text().toDouble();
    step_size = ui->plotMaxNStepSizeEdit->text().toDouble();
    int seed = ui->plotSeedEdit->text().toInt();

    if(step_size > min_deadline)
        deadline = step_size;
    else
        deadline = min_deadline;

    //if(step_size > min_deadline)
        //min_deadline = step_size;


    //output csv ---------------------------------------------------------------------------------------
    std::ofstream myfile;
    QString message;
    myfile.open (save_location.toUtf8().constData(), std::ios::trunc);
    if(myfile.is_open() == false)
    {
        ui->logEdit->insertPlainText(QString("cannot open output log file!\n"));
        setBusyState(false);
        return;
    }


    message = QString("Plot EMAC protocol, from d=%1 to d=%2, stepsize %3\n").arg(min_deadline)
            .arg(max_deadline).arg(step_size);
    myfile << message.toStdString().data();
    ui->logEdit->appendPlainText(message);
    myfile << "deadline[ms]; Numbers of nodes; calculation time [ms]; periods [ms]\n";
    myfile.close();

    for(; deadline <= max_deadline; deadline = deadline + step_size, k++)
    {
        double deadline_normalized = deadline/transmission_duration;

        //start timer
        timer.restart();

        for(n = seed; n < BUFFER_SIZE; n++)
        {
            if(EMAC_calculatePeriods((unsigned int) deadline_normalized,n,false) == false)
            {
                n--;
                time_elapsed = timer.elapsed();
                plot_storage.elapsed_time[k] = time_elapsed;

                if(EMAC_calculatePeriods((unsigned int) deadline_normalized,n,false) == false)  //failed twice -> cancel
                {
                    ui->logEdit->insertPlainText(QString("deadline = %1 ms  max n = %2\n").arg(deadline).arg(0));
                    plot_storage.number_of_nodes[k] = 0;
                    plot_storage.deadline[k] = (unsigned int) deadline;
                }
                else
                {
                    ui->logEdit->insertPlainText(QString("deadline = %1 ms  max n = %2\n").arg(deadline).arg(n));
                    plot_storage.number_of_nodes[k] = n;
                    plot_storage.deadline[k] = (unsigned int) deadline;
                }

                //only save increments of n
                if(n > nmax)
                {
                    nmax = n;
                    comp_time[comp_time_index] = time_elapsed;
                    comp_time_n[comp_time_index] = nmax;
                    comp_time_index++;
                }


                //copy period data
                for (int i = 0; i < shortest_selected_sequence_list.size(); ++i)
                {
                    //deadline_last = shortest_selected_sequence_list.at(i).deadline;
                    period_buf[i] = (unsigned int) shortest_selected_sequence_list.at(i).numbers[0];  // * ui->packetSizeEdit->text().toDouble();
                }

                //save in file (file has to be closed an reopened again to save changes immediately)
                myfile.open (save_location.toUtf8().constData(), std::fstream::out | std::fstream::app);
                if(myfile.is_open() == false)
                {
                    ui->logEdit->insertPlainText(QString("cannot open output log file!\n"));
                    goto fsw_end;
                }
                //message = QString("%1;%2\n").arg(plot_storage.deadline[k]).arg(plot_storage.number_of_nodes[k]);
                //myfile.write(message.toLocal8Bit(), message.count());
                myfile << plot_storage.deadline[k] << ";" << plot_storage.number_of_nodes[k] << ";" <<  time_elapsed << ";";
                //append periods
                for(unsigned int i = 0; i<n; i++)
                {
                    myfile << period_buf[i] * transmission_duration << ";";
                }
                myfile << "\n";
                myfile.close();
                myfile.close();

                //update GUI to show progress
                QCoreApplication::processEvents(QEventLoop::AllEvents, 10);
                QTextCursor c = ui->logEdit->textCursor();
                c.movePosition(QTextCursor::End);
                ui->logEdit->setTextCursor(c);

                break;
            }
        }
        seed = n;
    }

    plot_storage.step_size = step_size;
    plot_storage.steps = k;
    plot_storage.last_mode = 4; //EMAC

    ui->logEdit->insertPlainText(QString("finished!\n"));

    //save comp time in csv
    myfile.open (save_location.toUtf8().constData(), std::fstream::out | std::fstream::app);
    myfile << "\n\n" << "n; computation time;\n";
    for(int i=0; i<comp_time_index; i++)
    {
        myfile << comp_time_n[i] << ";" << comp_time[i] << ";\n";
    }
    myfile.close();

fsw_end:
    setBusyState(false);

    //plotSaveFromBufToFileButtonPressed();
}



void MainWindow::plotSaveFromBufToFileButtonPressed()
{
    //output csv ---------------------------------------------------------------------------------------
    std::ofstream myfile;
    QString message;
    myfile.open (save_location.toUtf8().constData(), std::ios::trunc);
    if(myfile.is_open() == false)
    {
        ui->logEdit->insertPlainText(QString("cannot open output log file!\n"));
        return;
    }

    message = QString("Plot ");

    if(plot_storage.last_mode == 3)
        message = QString("DEEP Analytic ");
    else if (plot_storage.last_mode == 1)
        message += "modified ";
    else if (plot_storage.last_mode == 2)
        message += "old ";
    else if (plot_storage.last_mode == 4)
       message += "EMAC ";
    else
        message += "unknown";

    message += QString("protocol, from d=%1 to d=%2, stepsize %3\n").arg(plot_storage.deadline[0])
            .arg(plot_storage.deadline[plot_storage.steps]).arg(plot_storage.step_size);
    myfile << message.toStdString().data();
    //save data in csv file
    myfile << "deadline[ms]; Numbers of nodes; calculation time [ms]\n";

    for(unsigned j = 0; j < plot_storage.steps; j++)
    {
        myfile << plot_storage.deadline[j] << ";" << plot_storage.number_of_nodes[j] << ";" << plot_storage.elapsed_time[j] <<  "\n";
    }
    //end output csv -------------------------------------------------------------------------------------

    ui->logEdit->insertPlainText(QString("saved succesfully in file\n"));
}



//////////////////
/// create period list functions
//////


void MainWindow::createPeriodListDEEPHeuristicButtonPressed()
{
    double nmin,nmax;
    bool seed_used;
    nmin = ui->createPeriodListNMinEdit->text().toDouble(); if(nmin<2) nmin=2;
    nmax = ui->createPeriodListNMaxEdit->text().toDouble();
    seed_used = ui->createPeriodListBox->isChecked();
    QString message;
    double deadline_last = 0;
    double deadline_last_normalized=0;
    double deadline_step_size=1; //in ms
    double time_base = 8 / ui->speedEdit->text().toDouble();
    double speed = ui->speedEdit->text().toDouble();

    //deadline list from previous sim
    double deadline_seed_list[30]=
        {1,
         3,
         5,
         11,
         17,
         27,
         38,
         47,
         66,
         85,
         110,
         120,
         180,
         190,
         220,
         280,
         350,
         380,
         450,
         490,
         530,
         590,
         710,
         740,
         790,
         930,
         990,
         1130,
         1290,
    };

    setBusyState(true);

    //read inputs
    updatePacketDeadlineList();

    //output csv ---------------------------------------------------------------------------------------
    std::ofstream myfile;
    myfile.open (save_location.toUtf8().constData(), std::ios::trunc);
    if(myfile.is_open() == false)
    {
        ui->logEdit->insertPlainText(QString("cannot open output log file!\n"));
        setBusyState(false);
        return;
    }

    message += QString("createPeriodListDEEPHeuristicButtonPressed n_min=%1  n_max=%2, seed_used=%3, speed=%4, lmax=%5 [bytes]\n").arg(nmin)
            .arg(nmax).arg(seed_used).arg(ui->speedEdit->text().toDouble()).arg(ui->packetSizeEdit->text().toDouble());

    ui->logEdit->insertPlainText(message);

    myfile << message.toStdString().data();
    myfile << "n; deadline; ; periods;\n";
    myfile.close();


    for(double n_iterator = nmin; n_iterator <= nmax; n_iterator++)
    {
        ui->logEdit->insertPlainText(QString("n=%1\n").arg(n_iterator));
        //debug
        //ui->logEdit->insertPlainText(QString("calculate for n=%1\n").arg(n_iterator));

        //calculate periods
        calculatePacketNumbers(n_iterator);//fill packet list

        if(seed_used == true)
        {
            if(n_iterator <= (sizeof(deadline_seed_list)/sizeof(double)))
            {
                //fill deadline list
                for(int i=0; i<n_iterator; i++)
                {
                    deadline_list_normalized[i] = deadline_seed_list[(int)n_iterator-2] /(8 / ui->speedEdit->text().toDouble());
                }
                deadline_last = deadline_seed_list[(int)n_iterator-2];
            }
            else
            {
                //slowly increase deadline until one was found
                while(1)
                {
                    //increase deadline
                    deadline_last += 10; //todo find a good step size

                    //fill buffers
                    for(int i=0; i<n_iterator; i++)
                    {
                        deadline_list_normalized[i] = deadline_last /(8 / ui->speedEdit->text().toDouble());
                    }

                    if(calculatePeriodsDEEPH(n_iterator) == true)
                    {
                        ui->logEdit->insertPlainText(QString("n=%1   deadline found: %2\n").arg(n_iterator).arg(deadline_last));
                        goto jumpi;
                    }
                }
            }
        }
        else
        {
            ui->logEdit->insertPlainText("not yet implemented -> activate seed check box\n");

            //TODO UNTESTED, check the following code
            for(double deadline_iterator=0; ; deadline_iterator += deadline_step_size)
            {
                //fill deadline list
                deadline_last = deadline_iterator;
                deadline_last_normalized = deadline_last / (8 / ui->speedEdit->text().toDouble());
                for(double i=0; i<n_iterator; i++)
                {
                    deadline_list_normalized[(int)i] = deadline_last /(8 / ui->speedEdit->text().toDouble());
                }

                //check and calculate periods
                if(calculatePeriodsDEEPH(n_iterator) == true)
                    break;
            }
        }




        calculatePeriodsDEEPH(n_iterator);

jumpi:
        //output results
        myfile.open (save_location.toUtf8().constData(), std::fstream::out | std::fstream::app);
        if(myfile.is_open() == false)
        {
            ui->logEdit->insertPlainText(QString("cannot open output log file!\n"));
            setBusyState(false);
            return;
        }
        myfile << n_iterator << ";" << deadline_last << ";" << ";";
        for(int i = 0; i<n_iterator; i++)
        {
            myfile << period_buf[i] << ";";
        }
        myfile << "\n";
        myfile.close();

        //update GUI to show progress
        QCoreApplication::processEvents(QEventLoop::AllEvents, 100);
        QTextCursor c = ui->logEdit->textCursor();
        c.movePosition(QTextCursor::End);
        ui->logEdit->setTextCursor(c);
    }

    setBusyState(false);
}


void MainWindow::createPeriodListDEEPAnalyticButtonPressed()
{
    double nmin,nmax;
    bool seed_used;
    nmin = ui->createPeriodListNMinEdit->text().toDouble();
    nmax = ui->createPeriodListNMaxEdit->text().toDouble();
    seed_used = ui->createPeriodListBox->isChecked();
    QString message;

    setBusyState(true);

    //read inputs
    updatePacketDeadlineList();


    //output csv ---------------------------------------------------------------------------------------
    std::ofstream myfile;
    myfile.open (save_location.toUtf8().constData(), std::ios::trunc);
    if(myfile.is_open() == false)
    {
        ui->logEdit->insertPlainText(QString("cannot open output log file!\n"));
        setBusyState(false);
        return;
    }

    message += QString("createPeriodListDEEPAnalyticButtonPressed n_min=%1  n_max=%2, seed_used=%3, speed=%4, lmax=%5 [bytes]\n").arg(nmin)
            .arg(nmax).arg(seed_used).arg(ui->speedEdit->text().toDouble()).arg(ui->packetSizeEdit->text().toDouble());

    ui->logEdit->insertPlainText(message);

    myfile << message.toStdString().data();
    myfile << "n; periods;\n";
    myfile.close();




    for(double n_iterator = nmin; n_iterator <= nmax; n_iterator++)
    {
        //debug
        //ui->logEdit->insertPlainText(QString("calculate for n=%1\n").arg(n_iterator));

        //calculate periods
        calculatePeriodsDEEPA(n_iterator);

        //output results
        myfile.open (save_location.toUtf8().constData(), std::fstream::out | std::fstream::app);
        if(myfile.is_open() == false)
        {
            ui->logEdit->insertPlainText(QString("cannot open output log file!\n"));
            setBusyState(false);
            return;
        }
        myfile << n_iterator << ";";
        for(int i = 0; i<n_iterator; i++)
        {
            myfile << period_buf[i] << ";";
        }
        myfile << "\n";
        myfile.close();
    }

    setBusyState(false);
}

void MainWindow::createPeriodListEMACButtonPressed()
{
    double nmin,nmax;
    bool seed_used;
    nmin = ui->createPeriodListNMinEdit->text().toDouble(); if(nmin<2) nmin=2;
    nmax = ui->createPeriodListNMaxEdit->text().toDouble();
    seed_used = ui->createPeriodListBox->isChecked();
    QString message;
    double deadline_last = 0;
    //double deadline_last_normalized=0;
    //double deadline_step_size=1; //in ms
    //double time_base = 8 / ui->speedEdit->text().toDouble();
    //double speed = ui->speedEdit->text().toDouble();

    //deadline list from previous sim
    double deadline_seed_list[30]=
    {1,
     1,
     4,
     6,
     12,
     17,
     30,
     40,
     50,
     60,
     80,
     90,
     100,
     130,
     150,
     170,
     190,
     220,
     250,
     280,
     300,
     330,
     400,
     430,
     460,
     500,
     540,
     600,
     660,
};

    setBusyState(true);

    //read inputs
    updatePacketDeadlineList();

    //output csv ---------------------------------------------------------------------------------------
    std::ofstream myfile;
    myfile.open (save_location.toUtf8().constData(), std::ios::trunc);
    if(myfile.is_open() == false)
    {
        ui->logEdit->insertPlainText(QString("cannot open output log file!\n"));
        setBusyState(false);
        return;
    }

    message += QString("createPeriodListEMACButtonPressed n_min=%1  n_max=%2, seed_used=%3, speed=%4, lmax=%5 [bytes]\n").arg(nmin)
            .arg(nmax).arg(seed_used).arg(ui->speedEdit->text().toDouble()).arg(ui->packetSizeEdit->text().toDouble());

    ui->logEdit->insertPlainText(message);

    myfile << message.toStdString().data();
    myfile << "n; deadline; ; periods;\n";
    myfile.close();


    for(double n_iterator = nmin; n_iterator <= nmax; n_iterator++)
    {
        //debug
        //ui->logEdit->insertPlainText(QString("calculate for n=%1\n").arg(n_iterator));

        //calculate periods
        calculatePacketNumbers(n_iterator);//fill packet list
        if(seed_used == true)
        {
            //fill deadline list
            for(int i=0; i<n_iterator; i++)
            {
                deadline_list_normalized[i] = deadline_seed_list[(int)n_iterator-2] /(8 / ui->speedEdit->text().toDouble());
            }
        }
        else
        {
            ui->logEdit->insertPlainText("not yet implemented -> activate seed check box\n");
        }




        EMAC_calculatePeriods((unsigned int) deadline_list_normalized[0], (unsigned int) n_iterator, false);

        //copy period data
        for (int i = 0; i < shortest_selected_sequence_list.size(); ++i)
        {
            deadline_last = shortest_selected_sequence_list.at(i).deadline;
            period_buf[i] = (unsigned int) shortest_selected_sequence_list.at(i).numbers[0];
        }

        //output results
        myfile.open (save_location.toUtf8().constData(), std::fstream::out | std::fstream::app);
        if(myfile.is_open() == false)
        {
            ui->logEdit->insertPlainText(QString("cannot open output log file!\n"));
            setBusyState(false);
            return;
        }
        myfile << n_iterator << ";" << deadline_last << ";" << ";";
        for(int i = 0; i<n_iterator; i++)
        {
            myfile << period_buf[i] << ";";
            //ui->logEdit->insertPlainText(QString("period_buf[%1]=%2\n").arg(i).arg(period_buf[i]));
        }
        myfile << "\n";
        myfile.close();

        //update GUI to show progress
        QCoreApplication::processEvents(QEventLoop::AllEvents, 100);
        QTextCursor c = ui->logEdit->textCursor();
        c.movePosition(QTextCursor::End);
        ui->logEdit->setTextCursor(c);
    }

    setBusyState(false);
}

//TODO
void MainWindow::plotMaxNPerDeadlinePKPressed() {

    //pkCalculatePeriods();
    //return;


    double min_deadline, max_deadline, step_size, deadline;
    int n, n_old = 2;
    unsigned int k=0;
    QElapsedTimer timer;
    long int time_elapsed;
    double comp_time[BUFFER_SIZE];
    int comp_time_n[BUFFER_SIZE];   //stores n for comp_time[]
    int comp_time_index = 0;

    setBusyState(true);

    std::ofstream myfile;
    QString message;
    myfile.open (save_location.toUtf8().constData(), std::ios::trunc);
    if(myfile.is_open() == false)
    {
        ui->logEdit->insertPlainText(QString("cannot open output log file!\n"));
        setBusyState(false);
        return;
    }
    min_deadline = ui->plotMaxNMinDEdit->text().toDouble();
    max_deadline = ui->plotMaxNMaxDEdit->text().toDouble();
    step_size = ui->plotMaxNStepSizeEdit->text().toDouble();

    if(step_size > min_deadline)
        deadline = 0;
    else
        deadline = min_deadline;


    message += QString("plotMaxNPerDeadlinePKPressed from d=%1 to d=%2, stepsize %3\n").arg(min_deadline)
            .arg(max_deadline).arg(step_size);

    message += "\n";
    myfile << message.toStdString().data();
    myfile << "deadline[ms]; Numbers of nodes; calculation time [ms]; periods\n";
    ui->logEdit->appendPlainText(message);

    //get information
    updatePacketDeadlineList();
    double transmission_duration = 8 / ui->speedEdit->text().toDouble();    //duration of time base (8bit)

    for(; deadline <= max_deadline; deadline = deadline + step_size, k++)
    {
        timer.restart();

        //update input for calculateMaxNPerD()
        for(unsigned int i = 0; i < BUFFER_SIZE; i++)
        {
            deadline_list_normalized[i] = (int) (deadline / transmission_duration);
        }

        n = pkCalculateMaxN(n_old, deadline, transmission_duration);

        time_elapsed = timer.elapsed();
        if(n > n_old)
        {
            n_old = n;
            comp_time[comp_time_index] = time_elapsed;
            comp_time_n[comp_time_index] = n_old;
            comp_time_index++;
            //std::cout << "n_highest: " << n_highest << "   time_elapsed: " << time_elapsed << std::endl;
        }
        ui->logEdit->insertPlainText(QString("deadline = %1 ms  max n = %2\n").arg(deadline).arg(n));
        myfile << deadline << ";" << n << ";" << time_elapsed <<  ";";
        //append periods
        for(int i=0; i<n; i++) {
            myfile << period_buf[i] << ";";
        }
        myfile  <<  "\n";

    }

    //save comp time in csv
    myfile << "\n\n" << "n; computation time;\n";
    for(int i=0; i<comp_time_index; i++)
    {
        myfile << comp_time_n[i] << ";" << comp_time[i] << ";\n";
    }
    myfile.close();
    ui->logEdit->appendPlainText("finished!");

    setBusyState(false);
}


int MainWindow::pkCalculateMaxN(int seed = 3, unsigned int max_deadline = 0, unsigned int packet_size_ = 0) {
    unsigned int n;
    unsigned int deadline;
    unsigned int packet_size;
    if(packet_size_ == 0)
        packet_size = ui->packetSizeEdit->text().toInt();
    else
        packet_size = packet_size_;
    double transmission_duration = 8 / ui->speedEdit->text().toDouble();    //duration of time base (8bit)

    //convert to normalized
    max_deadline = max_deadline / transmission_duration;

    for(n = seed; n < BUFFER_SIZE; n++)
    {

        //calculate periods for given n
        //calc used deadline with longest period
        //check if used deadline is greater than max deadline (if so, reduce n, if not repeat with n+1)


        deadline = pkMaxDeadline(n, packet_size);

        if(deadline > max_deadline) {
            n = n - 1;
            return n;
        }


        //update GUI to prevent GUI freezing
        QCoreApplication::processEvents(QEventLoop::AllEvents, 10);
    }
    return 0; //found n greater than BUFF size (should not happen)



/*  blockList = [];
    distances = [];
    t = 1;
    for(i = 1:N)
        k = 1;
        while(k<N)
            while(any(blockList == k*t))
                t = t + 1;
                k = 0;
            end;
            k = k + 1;
        end;
        distances(end+1) = t*(2*lmax);
        for(l = 1:N)
             blockList(end+1) = l*t;
        end;

    end;
*/
}

int MainWindow::pkMaxDeadline(int n, unsigned int packet_size_) {
    QList<int> blockList;
    QList<int> distances;
    int t = 1;
    int k;

    for(int i=0; i<n; i++) {
        k = 1;
        while(k<n) {
            while(blockList.indexOf(k*t) >= 0) {
                t++;
                k = 0;
            }
            k++;
        }
        distances.append(t*2*packet_size_);
        period_buf[i] = t*2*packet_size_;
        for(int j=1; j<=n; j++) {
            blockList.append(j*t);
        }
    }

    int longest_deadline = (n-1) * distances.at(n-1) + packet_size_;
    return longest_deadline;

    //debug: output data
    /*std::cout << "found periods n=" << n << std::endl;
    for(int i=0; i<n; i++) {
        std::cout << "" << distances.at(i) << " ; ";
    }
    std::cout << std::endl;
    std::cout << "longest deadline " << longest_deadline << std::endl;*/
}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Misc
///
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void MainWindow::clearText()
{
    ui->logEdit->clear();
}



void MainWindow::saveText()
{
    QString fileName = QFileDialog::getSaveFileName(this, tr("Save File"),
                                                 QDir::currentPath(),
                                                 tr("Text Files (*.txt)"));
    if (fileName.isEmpty()) return;

    QFile file(fileName);
    if (!file.open(QFile::WriteOnly | QFile::Text)) {
        ui->logEdit->appendPlainText(tr("Cannot write file %1:\n%2.")
                            .arg(fileName)
                            .arg(file.errorString()));

        return;
    }


    //save output
    QTextStream out(&file);
    out << ui->logEdit->toPlainText();
    file.close();
}

//called when pathEdit is changed
void MainWindow::updateSavePath()
{
    save_location = ui->pathEdit->text();

    //replace backslashes
    while(save_location.contains('\\')) {
        save_location.replace(save_location.indexOf('\\'), 1, '/');
    }
    ui->pathEdit->setText(save_location);
    //std::cout << "save_location: " << save_location.toLatin1().constData() << std::endl;
}

void MainWindow::aboutButtonPressed() {
    QMessageBox msgBox;
    msgBox.setText("<b>Inter-packet time finder for reliable communication</b>");
    msgBox.setInformativeText("Program to calculate inter-packet times (periods), max n per d and computation time"
" for the algorithms PK [1], DEEP Analytic [2], DEEP Heuristic [2] and eMAC [3][4]. Written by Philip Parsch (philip.parsch@gmail.com). Last updated"
" on 7 Jan 2022.\n\n"
"[1] under review\n"
"[2] On Enforcing Reliability in Unidirectional WSNs: A MAC-Based Approach (Parsch, Philip), PhD thesis, Department of Computer Science, TU Chemnitz, 2019.\n"
"[3] B. Andersson, N. Pereira, and E. Tovar. Delay-Bounded Medium Access for Unidirectional Wireless Links. In Proceedings of the International Conference "
"on Real-Time Networks and Systems (RTNS), 2007.\n"
"[4]  B. Andersson, N. Pereira, and E. Tovar. Delay-Bounded Medium Access for Unidirectional Wireless Links. Technical report, CISTER - Research Centre "
"in Real-Time and Embedded Computing Systems, 2007.\n");
    msgBox.exec();
}

void MainWindow::setBusyState(bool state) {
    if(state){
        //change label color to indicate busy state
        ui->busyLabel->setStyleSheet("QLabel { background-color : red; color : red; }");
    }
    else {
        //change label color to indicate idle state
        ui->busyLabel->setStyleSheet("QLabel { background-color : green; color : green; }");
    }
    //needed for gui to refresh
    QCoreApplication::processEvents(QEventLoop::AllEvents, 10);
}

