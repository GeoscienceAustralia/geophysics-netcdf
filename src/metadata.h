/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#ifndef _metadata_H
#define _metadata_H

#include "general_utils.h"
#include "file_utils.h"
#include "blocklanguage.h"

class cMetaDataTable{

public:
	std::vector<std::string> header;
	std::vector<std::vector<std::string>> records;

	cMetaDataTable(){ }

	cMetaDataTable(const std::string metadatafile, size_t nskip){

		FILE* fp = fileopen(metadatafile, "r");

		std::string str;
		std::vector<std::string> tokens;

		//Read header line
		filegetline(fp, str);
		split(str, '\t', tokens);
		header = tokens;
		size_t nfields = header.size();

		//Skip lines after header
		for (size_t i = 0; i < nskip; i++){
			filegetline(fp, str);
		}

		//Read remaining lines
		size_t k = 1 + nskip;
		while (filegetline(fp, str)){
			k++;
			tokens.resize(0);
			split(str, '\t', tokens);
			if (tokens.size() == nfields - 1){
				tokens.push_back("");
			}

			if (tokens.size() != nfields){
				std::printf("Error: On line %zu of file %s\n", k, metadatafile.c_str());
				std::printf("Error: The number of header items (%zu) does not match the number of data items (%zu)\n", nfields, tokens.size());
			}

			for (size_t i = 0; i < tokens.size(); i++){
				tokens[i] = trim(tokens[i]);
				tokens[i] = stripquotes(tokens[i]);
			}
			records.push_back(tokens);
		}
		fclose(fp);
	}

	cMetaDataTable(const cBlock& b, const std::string& id){
		std::vector<std::vector<std::string>> g = b.getblockleftright(id);
		for (size_t i = 0; i < g.size(); i++){
			addfield(g[i][0]);
			setfield(g[i][0], g[i][1]);
		}
	}

	bool addfield(const std::string fname){
		header.push_back(fname);
		for (size_t i = 0; i < records.size(); i++){
			const std::string empty;
			records[i].push_back(empty);
		}
		return true;
	}

	bool setfield(const std::string fname, const std::string value, const size_t recindex){
		int k = findkeyindex(fname);
		if (k < 0)return false;
		records[recindex][(size_t)k] = value;
		return true;
	}

	bool setfield(const std::string fname, const std::string value){
		int k = findkeyindex(fname);
		if (k < 0)return false;
		if (records.size() == 0){
			records.resize(1);
			records[0].resize(header.size());
		}

		for (size_t i = 0; i < records.size(); i++){
			records[i][(size_t)k] = value;
		}
		return true;
	}

	int findkeyindex(const std::string& fname){
		for (size_t i = 0; i < header.size(); i++){
			if (strcasecmp(header[i].c_str(),fname.c_str())==0){
				return (int)i;
			}
		}
		return -1;
	}

	std::vector<size_t> findmatchingrecords(const size_t keyindex, const int value){
		std::vector<size_t> indices;
		for (size_t i = 0; i < records.size(); i++){
			int n = atoi(records[i][keyindex].c_str());
			if (n == value){
				indices.push_back(i);
			}
		}
		return indices;
	}

	std::vector<size_t> findmatchingrecords(const std::string key, const int value)	{
		int keyindex = findkeyindex(key);
		if (keyindex < 0){
			std::vector<size_t> indices;
			return indices;
		}
		return findmatchingrecords((size_t)keyindex, value);
	}

	void printrecord(const size_t n){
		for (size_t mi = 0; mi < header.size(); mi++){
			std::printf("%s: %s\n", header[mi].c_str(), records[n][mi].c_str());
		}
		std::printf("\n");
	}
};

class cMetaDataRecord{

public:
	std::vector<std::string> header;
	std::vector<std::string> values;

	cMetaDataRecord(){};

	cMetaDataRecord(const cBlock& b, const std::string& id){
		cMetaDataTable T(b, id);
		header = T.header;
		values = T.records[0];
	}

	int findkeyindex(const std::string& fname){
		for (size_t i = 0; i < header.size(); i++){
			if (header[i] == fname){
				return (int)i;
			}
		}
		return -1;
	}

	void print(){
		for (size_t mi = 0; mi < header.size(); mi++){
			std::printf("%s: %s\n", header[mi].c_str(), values[mi].c_str());
		}
		std::printf("\n");
	}
};

#endif

