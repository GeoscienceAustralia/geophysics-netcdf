/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2016.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#include <cstdio>
#include <vector>
#include <limits>

#include <netcdf>
#include "general_utils.h"
#include "file_utils.h"
#include "geophysics_netcdf.h"

NcType nctype(const short){ return ncShort; }
NcType nctype(const int){ return ncInt; }
NcType nctype(const unsigned int){ return ncUint; }
NcType nctype(const float){ return ncFloat; }
NcType nctype(const double){ return ncDouble; }
NcType nctype(const std::string){ return ncString; }

size_t cGeophysicsVar::line_index_start(const size_t& index){
	size_t start = get_parent()->get_line_index_start(index);	
	return start;
}

size_t cGeophysicsVar::line_index_count(const size_t& index){
	size_t count = get_parent()->get_line_index_count(index);	
	return count;
}