#if defined(_WIN32)
	//#include"dirent.h"
	#include <conio.h>
#else
	#include <dirent.h>
	#include <math.h>
	#include <cstring>
	#include <stdlib.h>
#endif
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <complex>
#include <time.h>
#include "mpi.h"

#define DEFAULT_SHIFT 1e-6
#define EXIT_MESSAGE "\nPress any key to exit"
#define W_I (de*i_i)
#define PRINT_INTERVAL 10
#define WORK_TIME (time(NULL)-::startTime)

using std::vector;
using std::string;
using std::cout;
using std::cerr;
using std::endl;

typedef std::complex<double> dcomplex;

const double MyPI=2*asin(1.0);
const time_t startTime=time(NULL);

int proc_number;//номер процессора
int n;//количество процессоров
vector<unsigned int>steps_per_proc;//количество шагов в одном процессе
//int my_error_code=0;

template<class T>
inline string tostring(//Преобразование произвольного типа, для которого определена операция <<, в string
	const T&Tvalue
	)
{
	std::ostringstream ovalue;
	ovalue<<Tvalue;
	return ovalue.str();
}

void print_progress(string string_to_print)
{
	//Вывод на экран с интервалом
	time_t runTime=WORK_TIME;//time(NULL);
	static time_t previous_runTime=0;//Время работы, сохраненное с предыдущего вывода на экран
	static bool printIsPermitted=true;
	if(runTime-previous_runTime>=PRINT_INTERVAL)//Если время работы цикла с предыдущей итерации превышает время, заданное в PRINT_INTERVAL (в секундах), ...
	{
		printIsPermitted=true;//то разрешаем вывод на экран. Первый вывод всегда осуществляется
		previous_runTime=runTime; // Перезаписываем время работы, сохраненное с предыдущего вывода на экран
	}
	if(printIsPermitted)
	{
		cerr<<string_to_print;
		printIsPermitted=false;
	}
}

inline vector<double>vector_complex2double(const vector<dcomplex>&a)
{
	vector<double>r(a.size());
	for(unsigned int i=0;i<r.size();i++)
		r[i]=a[i].real();
	return r;
}

inline int Find_Proc_Number(//Функция, которая определяет принадлежность i_r какому-либо процессу
	int i_r,//(Z)
	int size,//(Z_length) количество i_r (элементов массива)
	int n//(NCELLS) количество процессов
	)
{
	int scale=size/n;
	int cell_number=i_r/scale;
	//if(cell_number<0)      cell_number=0;
	if(cell_number>=n)
		cell_number=n-1;
	return cell_number;//возвращаемое значение - искомый номер элементарной ячейки, нумерация от нуля
}//Эта функция не будет работать, если количество элементов массива меньше количества процессов.
//Исправить это!!

/*	Calculate the Kramers-Kronig transformation on imaginary part of dielectric

	Doesn't correct for any artefacts resulting from finite window function.

	Args:
		de (float): Energy grid size at which the imaginary dielectric constant
			is given. The grid is expected to be regularly spaced.
		eps_imag (np.array): A numpy array with dimensions (n, 3, 3), containing
			the imaginary part of the dielectric tensor.
		cshift (float, optional): The implemented method includes a small
			complex shift. A larger value causes a slight smoothing of the
			dielectric function.

	Returns:
		A numpy array with dimensions (n, 3, 3) containing the real part of the
		dielectric function.*/
	//Отсюда: https://github.com/utf/kramers-kronig ; https://github.com/utf/kramers-kronig/blob/master/kkr.ipynb
vector<double>kkr
	(
		double de,
		const vector<double>&eps_imag,
		double shift=DEFAULT_SHIFT,
		bool direction=true // true - Действительная часть -> мнимая (r2i)
							//false - Мнимая часть -> действительная (i2r)
	)
{
	vector<dcomplex>eps_real;
	eps_real.resize(eps_imag.size(),0.0);
	for (unsigned int i_r=0;i_r<eps_imag.size();i_r++)
	{
		if (Find_Proc_Number(i_r,eps_imag.size(),::n)< ::proc_number)
			continue;
		if (Find_Proc_Number(i_r,eps_imag.size(),::n)> ::proc_number)
			break;
		double w_r=de*i_r;
		double w_r2=w_r*w_r;
		dcomplex total(0.0,0.0);
		const dcomplex c_shift(0.0,shift);
		unsigned int i_i=1;
		if(direction)
		{
			for (;i_i<eps_imag.size()-1;i_i++) // Во внутреннем цикле условные операторы не используются, чтобы увеличить быстродействие
			{
				double w_i = W_I;
				total+= eps_imag[i_i]/(w_r2-w_i*w_i+c_shift);//w_r вынесена из внутреннего цикла во внешний для увеличения быстродействия
			}
			total+=(eps_imag[0]/(w_r2+c_shift)+eps_imag[i_i]/(w_r2-W_I*W_I+c_shift))/2.;
			//eps_real.push_back(2.0*total*w_r*de/MyPI);
		}
		else
		{
			for (;i_i<eps_imag.size()-1;i_i++)
			{
				double w_i = W_I;
				total+= eps_imag[i_i]*w_i/(w_r2-w_i*w_i+c_shift);
			}
			total+=eps_imag[i_i]*W_I/(w_r2-W_I*W_I+c_shift)/2.;
			//eps_real.push_back(2.0*total*(-1.0)*de/MyPI);
		}
		eps_real[i_r]=2.0*total*(direction?w_r:-1.0)*de/::MyPI;
		::steps_per_proc[::proc_number]++;// В соответствующий элемент массива заносим, что процесс посчитал один шаг

		if(!::proc_number)
			print_progress
			(
				tostring(WORK_TIME)+(string)"s: .."+tostring(100.0*i_r*::n/((double)eps_imag.size()))+(string)"% is done..\n"
			);
	}
	return vector_complex2double(eps_real);
}

void error_quit (string message_to_print="Unknown error\n", int error_code=-1)
{
	if(!::proc_number)
	{
		cout<<message_to_print;
		#if defined(_WIN32)
			cerr<<EXIT_MESSAGE<<endl;
			_getch();
		#endif
	}
	//MPI_Barrier(MPI_COMM_WORLD);
	//MPI_Finalize();
	exit(error_code);
	//::my_error_code=error_code;
	//for(unsigned int this_proc_number=1;this_proc_number<::n;this_proc_number++)
		//MPI_Send(&::my_error_code,1,MPI_INT,this_proc_number,0,MPI_COMM_WORLD);
}//Сделать корректное завершение всех процессов в случае ошибки!!
//Сейчас в случае ошибки завершается только первый процесс, остальные могут остаться, если не сработает система контроля очереди на кластере

int main(int argc, char * argv[])
{
	//cerr<<"%"<<endl;
	MPI_Init(&argc,&argv); //инициализация MPI
	MPI_Comm_rank(MPI_COMM_WORLD,&::proc_number); //Определение номера текущего процесса
	MPI_Comm_size(MPI_COMM_WORLD,&::n); //Определение количества процессоров
	if (!::proc_number)
		cerr<<WORK_TIME<<"s: "<<"Program is started, "<<::n<<" process(es) run(s).."<<endl;

	//vector<double>Im;//Скорее всего, этот массив не нужен. Результат можно сохранить в cond для экономии памяти
	vector<double>cond;// Сюда считываем значения из файла
	double de;
	double shift;
	bool direction=true;// true - Действительная часть -> мнимая (r2i)
						//false - Мнимая часть -> действительная (i2r)
	string firststring;
	double e0;
	unsigned int cond_size;
	if (!::proc_number)
	{
		if((argc<2)||(argc>5))
			error_quit((string)"Wrong arguments!\nThe syntax is:\n"+(string)(argv[0])+(string)" <File> [<shift(default="+tostring(DEFAULT_SHIFT)+(string)")>] [<r2i(default) or i2r>] [tail_length(default=0)in_grid_nodes]\n",
			1);

		std::ifstream infile(argv[1]);
		if(infile.fail())
			error_quit((string)"Couldn't open input file \""+tostring(argv[1])+(string)"\"!\n",
			2);
		cerr<<WORK_TIME<<"s: "<<"File \""<<argv[1]<<"\" is opened."<<endl;

		shift=argc>2?atof(argv[2]):DEFAULT_SHIFT;

		if(argc>3)
			if(!strcmp(argv[3],"r2i"))
				direction=true; // Направление преобразования: Действительная часть -> мнимая (r2i)
			else if(!strcmp(argv[3],"i2r"))
				direction=false;// Мнимая часть -> действительная (i2r)
			else
				error_quit((string)"Unknown parameter \""+tostring(argv[3])+(string)"\"!\n",
				3);

		bool FirstStringIsPassed=false;
		bool SecondStringIsPassed=false;
		bool ThirdStringIsPassed=true;
		double e1;
		while(!infile.eof())
		{
			string row;
			double first, second;
			getline(infile,row);
			std::istringstream(row)>>first>>second;
			if(!FirstStringIsPassed)
			{
				FirstStringIsPassed=true;
				firststring=row;
				continue;
			}
			if(!ThirdStringIsPassed)
			{
				e1=first;
				ThirdStringIsPassed=true;
			}
			if(!SecondStringIsPassed)
			{
				e0=first;
				SecondStringIsPassed=true;
				ThirdStringIsPassed=false;
			}
			if(row.length())
				cond.push_back(second); //Пропускаем строки нулевой длины. Такая строка обычно бывает в конце входного файла
		}
		infile.close();
		cerr<<WORK_TIME<<"s: "<<"File \""<<argv[1]<<"\" is read."<<endl;
		cerr<<"Array of conductivity values of "<<cond.size()<<" elements is created."<<endl;
		de=e1-e0;

		if (argc>4)
		{
			unsigned int tail_length=atoi(argv[4]);
			for(unsigned int i=0;i<tail_length;i++)
				cond.push_back(cond.back());
			cerr<<WORK_TIME<<"s: "<<"Array of conductivity values is enlarged with "<<tail_length<<" elements."<<endl;
		}
		cond_size=cond.size();

		cerr<<WORK_TIME<<"s: "<<"Kramers-Kronig transform is started..."<<endl;
	}
	MPI_Bcast(&de,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&cond_size,1,MPI_INT,0,MPI_COMM_WORLD);
	if(::proc_number)
		cond.resize(cond_size,0.0);
	MPI_Bcast(&cond[0],cond_size,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&shift,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&direction,1,MPI_INT,0,MPI_COMM_WORLD);

	::steps_per_proc.resize(::n,0);//количество шагов в одном процессе
	cond=kkr(de,cond,shift,direction);//Преобразование Крамерса-Кронига
	vector<double>cond_all(cond_size,0.0);
	vector<unsigned int>steps_per_proc_all(::n,0);
	MPI_Allreduce(&cond[0],&cond_all[0],cond_size,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(&::steps_per_proc[0],&steps_per_proc_all[0],::n,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
	
	if (!::proc_number)
	{
		cout<<firststring<<endl;
		for(unsigned int i=0;i<cond_all.size();i++)
			cout<<e0+de*i<<"\t"<<cond_all[i]<<endl;
		cerr<<"Total number of conductivity array elements: "<<cond_all.size()<<endl;
		cerr<<"Distribution i_r steps per processes:"<<endl;
		for (unsigned int this_proc_number=0;this_proc_number<steps_per_proc_all.size();this_proc_number++)
			cerr<<steps_per_proc_all[this_proc_number]<<' ';
		cerr<<endl;
		cerr<<WORK_TIME<<"s: "<<"All done!!"<<endl;
	}
	
	MPI_Finalize();
	return 0;
}
