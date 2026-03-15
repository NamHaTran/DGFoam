#ifndef DGCONTROLLER_H_INCLUDED
#define DGCONTROLLER_H_INCLUDED
#include <string>

/*Main controller function*/
void Executer();

void checkCommandLine(std::string cmd);

/*Processing*/
void Processing();

/*PreProcessing*/
void PreProcessing();

/*PostProcessing*/
void PostProcessing();

#endif // DGCONTROLLER_H_INCLUDED
