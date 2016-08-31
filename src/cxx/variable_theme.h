/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#ifndef _variable_theme_H
#define _variable_theme_H

#include "blocklanguage.h"
#include "intrepid.h"

class cVariable{

	ILDataset* pDataset;

public:
	std::string standardname;
	std::string longname;
	std::string aliaskey;
	std::string aliasvalue;
	std::vector<std::string> likelyfieldnames;
	std::vector<std::string> fieldnames;
	std::vector<std::string> requiredbythemes;
	std::string units;

	cVariable(cBlock b){
		standardname = b.getstringvalue("StandardName");
		longname = b.getstringvalue("LongName");
		aliaskey = b.getstringvalue("SurveyInfoAlias");
		std::string s = b.getstringvalue("LikelyNames");
		likelyfieldnames = tokenize(s);
		units = b.getstringvalue("Units");
	}

	size_t setdataset(ILDataset* pDa){

		pDataset = pDa;
		fieldnames.resize(0);
		if (aliaskey.size() > 0 && aliaskey != cBlock::ud_string()){
			aliasvalue = pDataset->surveyinfofieldname(aliaskey);
			if (aliasvalue.size() > 0){
				fieldnames.push_back(aliasvalue);
			}
		}

		for (size_t i = 0; i < likelyfieldnames.size(); i++){
			std::string s = pDataset->fieldnamelike(likelyfieldnames[i]);
			if (s.size()>0){
				bool alreadyhavefield = false;
				for (size_t j = 0; j<fieldnames.size(); j++){
					if (s == fieldnames[j]){
						alreadyhavefield = true;
					}
					break;
				}
				if (alreadyhavefield == false){
					fieldnames.push_back(s);
				}
			}
		}
		return fieldnames.size();
	}

	bool hasaliaskey(){
		if (aliaskey == cBlock::ud_string()){
			return false;
		}
		return true;
	}

	bool aliasfieldexists(){
		if (hasaliaskey()){
			std::string s = pDataset->fieldnamelike(aliasvalue);
			if (s.size() > 0){
				return true;
			}
		}
		return false;
	}

	bool reportproblems(){

		bool status = true;
		if (hasaliaskey()){
			if (pDataset->hassurveyinfoid(aliaskey)){
				if (aliasfieldexists() == false){
					std::printf("Warning 3: %s\n", standardname.c_str());
					std::printf("\tAlias key %s is defined but field %s does not exist\n", aliaskey.c_str(), aliasvalue.c_str());
					status = false;
				}
			}
		}

		if (fieldnames.size() > 1){
			std::printf("Warning 4: %s\n", standardname.c_str());
			std::printf("\tMore than one likely field was found: ");
			for (size_t i = 0; i < fieldnames.size(); i++){
				std::printf("%s ", fieldnames[i].c_str());
			}
			std::printf("\n");
			status = false;
		}
		return status;

	}

	std::string likelyfield(){
		return fieldnames[0];
	}
};

class cTheme{

public:
	bool process;
	std::string Name;
	std::vector<std::vector<std::string>> required;

	cTheme(cBlock b){
		Name = b.getstringvalue("Name");
		process = b.getboolvalue("Process");
		if (process == false)return;

		std::vector<std::string> v = b.getblockstrings("RequiredVariables");
		for (size_t i = 0; i < v.size(); i++){
			required.push_back(tokenize(v[i]));
		}
	};

};

#endif

