#include <gflags/gflags.h>
#include <iostream>
#include <time.h>
#include <string>
#include <stdlib.h>
#include <sstream>

DEFINE_string(tex_target, "None", "LaTeX target file");
DEFINE_string(message, "None", "Sets commit message");
DEFINE_string(expt, "None", "Sets current experiment");

int main(int argc, char *argv[]){
	using namespace std;
	google::ParseCommandLineFlags(&argc, &argv, true);

	string target (FLAGS_tex_target);
	string message (FLAGS_message);
	string expt (FLAGS_expt);

	string key = "None";

	long timeinfo = time(0);
	stringstream t;
	t << timeinfo;
	string currtime (t.str());

	if (target.compare(key) != 0){
		string fmt;
		fmt = "pdflatex " + target;
		system(fmt.c_str());
	}
	else{;}

	if (expt.compare(key) != 0){
		string fmt;
		fmt = "git add " + expt;
		system(fmt.c_str());
	}
	else{
		system("git add .");
	}

	string output;
	output = system("git status");
	cout << output << endl;
	cout << "continue to Continue" << endl;

	string confkey ("continue");
	string conf;
	cin >> conf;
	if (conf.compare(confkey) == 0){
		if (message.compare(key) != 0){
			string fmt;
			fmt = "git commit -m " + message;
			system(fmt.c_str());
		}
		else{
			string fmt;
			fmt = "git commit -m " + currtime;
			system(fmt.c_str());
		}
	}
	system("git push -u origin master");
}
