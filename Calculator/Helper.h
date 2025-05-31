#ifndef HELPER
#define HELPER

#include <string>
#include <vector>
class Helper
{
	std::vector<std::string> output;
	
public:

	Helper();
	~Helper();
	void const addCommand(std::string input);

};

#endif