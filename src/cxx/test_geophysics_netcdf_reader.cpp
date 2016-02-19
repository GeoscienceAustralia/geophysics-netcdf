/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#include <cstdio>
#include <netcdf>
#include <vector>
#include <limits>

#define USEGLOBALSTACKTRACE
#ifdef USEGLOBALSTACKTRACE
#include "stacktrace.h"
cStackTrace globalstacktrace;
#endif

using namespace netCDF;
using namespace netCDF::exceptions;
#include "geophysicsncfile.h"
#include "general_utils.h"
#include "file_utils.h"
#include "blocklanguage.h"

int main(int argc, char** argv)
{
	_GSTITEM_
	//if (argc != 2){
	//	printf("Usage: %s nc_file_path\n", argv[0]);
	//	return 1;
	//}

	try
	{
		//std::string ncpath = argv[1];
		std::string ncpath = "http://dapds00.nci.org.au/thredds/dodsC/uc0/rr2_dev/rcb547/magrad_tests_indexed_v2/GSQP1029MAG.nc";
		//std::string ncpath   = "http://dapds00.nci.org.au/thredds/dodsC/uc0/rr2_dev/rcb547/magrad_tests_indexed_v2/P583MAG.nc";
		

		cGeophysicsNcFile ncfile(ncpath, NcFile::read);
		
		std::vector<int> linenumber;
		ncfile.getLineNumbers(linenumber);

		std::vector<int> flightnumber;
		ncfile.getFlightNumbers(flightnumber);

		std::vector<double> v1,v2,v3;
		for (size_t i = 0; i < ncfile.nlines(); i=i+100){			
			double t1 = gettime();			
			ncfile.getVarByLineIndex("mag_microLevelled", v1, i);
			ncfile.getVarByLineIndex("latitude_GDA94", v2, i);
			ncfile.getVarByLineIndex("longitude_GDA94", v3, i);
			double t2 = gettime();
			printf("%lu nsamples=%lu time=%lf\n", i, v1.size(), t2 - t1);
		}
		prompttocontinue();
		
	}

	catch (NcException& e)
	{
		_GSTPRINT_
		std::cout << e.what() << std::endl;
		return 1;
	}
	catch (std::exception& e)
	{
		_GSTPRINT_
		std::cout << e.what() << std::endl;
		return 1;
	}
	std::printf("Success\n");
	return 0;  // successfully terminated
}


