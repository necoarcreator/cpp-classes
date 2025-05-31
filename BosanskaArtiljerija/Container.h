#ifndef CONT
#define CONT

#define DEBUG

#include <string>
#include <list>
#include <tuple>
#include <utility>
//#include <vector>

using namespace std;
struct Container //этот класс хранит в себе все возможные данных, соотносит вводные с ними
{
	//потом надо добавить минимальные расстояния
	list<tuple<string, double, size_t, double, double>> artyRangesAmmo; //название артустановок, дальность стрельбы, число снарядов в б/к, скорострельность заявленная, калибр 
	list<string> shellTypes; //возможные типы зарядов боеприпасов
	list<string> shellFuzes;  //возможные виды взрывателей

	list<pair<string, size_t>> targetsPriorities; //список возможных целей и их приоритеты
	list<string> missions; //возможные огневые задачи
	list<string> targetGeometries; //возможные типы целей по геометрии
	list<string> preparingTypes; //возможные типы пристрелки
	
		//
	//vector<pair<double, double>> artyCoords;

	Container(string file1, string file2);
#ifdef DEBUG
	Container(); //для теста
#endif
	~Container();

	bool checkCorrectArty(string _artyType, size_t _numInBattery, size_t _numAmmo, string _shellType,
		string _shellFuzes, double _temperature); //проверить, есть ли такой вид артиллерии
	bool checkCorrectTarget(string _targetName, string _mission, string _entrenchedArmoured, string _targetGeom, string _preparingType); //проверить, есть ли такой вид цели
	size_t getPriority(string _targetName);
	//double measureFireRate(size_t _howManyShells, double _temperature, string _artyName, string _shellType); //рассчитать фактическую скорострельность

	size_t getNormStationary(pair<string, size_t> _caliber, bool _isCluster = false, 
		bool _isDoubleIndex = false, bool _isShortened = false, bool _isEntrenchedOrArmored = false,  
		bool _isDestruction = false, bool _isActiveReactive = false, double range = 5.0, size_t _numInBattery = 1, string _targetName = "InfOpened");

	size_t getNormMoving(pair<string, size_t> _caliber, bool _isCluster, bool _isArmoured, size_t _numInBattery, string _targetName = "column");
	size_t getNormFireAtWill(pair<string, size_t> _caliber, double _meters, double _minutes);
	size_t getNormSmoke(pair<string, size_t> _caliber, double _meters, double _minutes, bool _isSideWind);

	double getRange(string _artyType);
	
	double getCaliber(string _artyType);

	string getType(string _artyType);

	bool checkStationary(string _targetType);

	size_t measureFireRate(pair<string, size_t> _caliber, size_t _charge, double _temperature, size_t _numShells, string _specificType);

	string getSpecificType(string _artyType);

	size_t getCharge(string _artyName, string _targetName, string _fuzeType, bool _isEntrenchedOrArmored, double _range);
};


#endif