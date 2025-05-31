#include "Container.h"

#define DEBUG

#include <string>
#include <list>
#include <tuple>
#include <vector>
#include <map>
#include <iostream>
#include <algorithm>
#include <memory>
#include <limits.h>
#include <exception>

using namespace std;
using std::get;

#ifdef DEBUG

Container::Container()
{
	artyRangesAmmo = { {"2C1", 15.2, 40, 4, 122.0}, {"2C3", 20.5, 46, 3, 152.0}, {"2C12B", 7.1, 60, 9, 82.0}};

	shellTypes = {"HE", "shrapnel", "cluster"};

	shellFuzes = {"percussion", "radio", "remoteControl"};

	targetsPriorities = { {"CP",1}, {"tank", 2}, {"infEntranched", 3}, {"infOpened", 4} };

	missions = { "destruction", "supression", "harassment" };

	targetGeometries = { "singular", "column", "squared" };

	preparingTypes = { "full", "eye", "shortened" };
}

#endif

Container::Container(string file1, string file2) { cout << "i'm dummy" << endl; }
Container::~Container() { cout << "that's all folks!" << endl; }

bool Container::checkCorrectArty(string _artyType, size_t _numInBattery, size_t _numAmmo, string _shellType,
	string _shellFuzes, double _temperature)
{
	  
	auto it = find_if(artyRangesAmmo.begin(), artyRangesAmmo.end(), [&](auto const& _n)
		{
			return get<0>(_n) == _artyType; //если артсистема есть в списке

		});

	 bool first = (it != artyRangesAmmo.end()) ? true : false;
	 size_t reducedNumShells = static_cast<size_t>(_numAmmo / _numInBattery);
	 bool second = ((reducedNumShells > 0)&&(reducedNumShells < get<2>(*it))) ? true : false; //если снарядов больше 0 но меньше максимума
	 bool third = (find(shellTypes.begin(), shellTypes.end(), _shellType) != shellTypes.end()) ? true : false;
	 bool fourth = (find(shellFuzes.begin(), shellFuzes.end(), _shellFuzes) != shellFuzes.end()) ? true : false;
	 bool fifth = (_temperature > -50.0 && _temperature < +50.0) ? true : false;
	 
	 try
	 {
		 if (first && second && third && fourth && fifth) { return true; } //если запрос корректный
		 else
		 {
			 throw 2;
		 }
	 }
	 catch (int e)
	 {
		 cerr << "A problem with info in artillery list occured. Try check the file." << endl;
		 exit(e);
	 }
}
bool Container::checkCorrectTarget(string _targetName, string _mission, string _entrenchedArmoured, string _targetGeom, string _preparingType)
{

	auto it = find_if(targetsPriorities.begin(), targetsPriorities.end(), [&](auto const& _n)
		{
			return get<0>(_n) == _targetName; //если артсистема есть в списке

		});

	bool first = (it != targetsPriorities.end()) ? true : false;

	bool second = (find(missions.begin(), missions.end(), _mission) != missions.end()) ? true : false; //если снарядов больше 0 но меньше максимума
	bool third = (_entrenchedArmoured == "-" || _entrenchedArmoured == "entrenched" || _entrenchedArmoured == "armoured") ? true : false;
	bool fourth = (find(targetGeometries.begin(), targetGeometries.end(), _targetGeom) != targetGeometries.end()) ? true : false;
	bool fifth = (find(preparingTypes.begin(), preparingTypes.end(), _preparingType) != preparingTypes.end()) ? true : false;

	try
	{
		if (first && second && third && fourth && fifth) { return true; } //если запрос корректный
		else
		{
			throw 2;
		}
	}
	catch (int e)
	{
		cerr << "A problem with info in target list occured. Try check the file." << endl;
		exit(e);
	}
}
size_t Container::getPriority(string _targetName)
{
	try
	{
		auto it = find_if(targetsPriorities.begin(), targetsPriorities.end(), [&](auto const& _n)
			{return _n.first == _targetName; });

		auto it2 = targetsPriorities.begin();

		//if (it != targetsPriorities.end())
		//{
			//advance(it2, distance(targetsPriorities.begin(), it));
			//return it->second; //вернуть приоритет такой-то цели
		//}
		//else { throw 1; }
	}
	catch (int e)
	{
		cerr << "Error while searching priority for a target!\n";
		exit(e);
	}

}

size_t Container::getNormStationary(pair<string, size_t> _caliber, bool _isCluster,
	bool _isDoubleIndex, bool _isShortened, bool _isEntrenchedOrArmored, bool _isDestruction, bool _isActiveReactive, double _range, size_t _numInBattery, string _targetName)
{
	map<vector<string>, size_t> targets = { {{"artyEntrenched", "mortarEntranched"}, 0},
											{{"radarsAuto", "radioAuto", "AAOpened"}, 0 },
											{{"infEntranched", "CP", "tank", "APC"}, 0},
											{{"infOpened"}, 0},
											{{"autoCPOpened"}, 0},
											{{"ATGMOpened", "ATGunOpened", "otherUnArmouredOpened"}, 0} };
	

	auto it = targets.begin();
	if (_caliber.first == "howitzer")
	{
		switch (_caliber.second)
		{
		case 76:
			it->second = 540;
			advance(it, 1);
			it->second = 450;
			advance(it, 1);
			it->second = 450;
			advance(it, 1);
			it->second = 90;
			advance(it, 1);
			it->second = 150;
			advance(it, 1);
			it->second = 900;
			break;
		case 85:
			it->second = 480;
			advance(it, 1);
			it->second = 400;
			advance(it, 1);
			it->second = 400;
			advance(it, 1);
			it->second = 85;
			advance(it, 1);
			it->second = 120;
			advance(it, 1);
			it->second = 800;
			break;
		case 100:
			it->second = 360;
			advance(it, 1);
			it->second = 300;
			advance(it, 1);
			it->second = 300;
			advance(it, 1);
			it->second = 55;
			advance(it, 1);
			it->second = 80;
			advance(it, 1);
			it->second = 350;
			break;
		case 120:
			it->second = 200;
			advance(it, 1);
			it->second = 160;
			advance(it, 1);
			it->second = 150;
			advance(it, 1);
			it->second = 35;
			advance(it, 1);
			it->second = 50;
			advance(it, 1);
			it->second = 300;
			break;
		case 122:
			it->second = 240;
			advance(it, 1);
			it->second = 200;
			advance(it, 1);
			it->second = 180;
			advance(it, 1);
			it->second = 40;
			advance(it, 1);
			it->second = 50;
			advance(it, 1);
			it->second = 300;
			break;
		case 130:
			it->second = 220;
			advance(it, 1);
			it->second = 180;
			advance(it, 1);
			it->second = 160;
			advance(it, 1);
			it->second = 40;
			advance(it, 1);
			it->second = 50;
			advance(it, 1);
			it->second = 300;
			break;
		case 152:
			if (_isCluster)
			{
				it->second = 60;
				advance(it, 1);
				it->second = 50;
				advance(it, 1);
				it->second = UINT_MAX;
				advance(it, 1);
				it->second = 8;
				advance(it, 1);
				it->second = 15;
				advance(it, 1);
				it->second = 100;
			}
			else
			{
				it->second = 180;
				advance(it, 1);
				it->second = 150;
				advance(it, 1);
				it->second = 120;
				advance(it, 1);
				it->second = 25;
				advance(it, 1);
				it->second = 40;
				advance(it, 1);
				it->second = 300;

			}
			break;
		case 203:
			if (_isCluster)
			{
				it->second = 60;
				advance(it, 1);
				it->second = 30;
				advance(it, 1);
				it->second = UINT_MAX;
				advance(it, 1);
				it->second = 7;
				advance(it, 1);
				it->second = 10;
			}
			else
			{
				it->second = 180;
				advance(it, 1);
				it->second = 100;
				advance(it, 1);
				it->second = 40;
				advance(it, 1);
				it->second = 20;
				advance(it, 1);
				it->second = 30;

			}
			advance(it, 1);
			it->second = UINT_MAX;
			break;
		}
	}
	else if (_caliber.first == "mortar")
	{
		switch (_caliber.second)
		{
		case 82:
			it->second = UINT_MAX;
			advance(it, 1);
			it->second = UINT_MAX;
			advance(it, 1);
			it->second = 700;
			advance(it, 1);
			it->second = 95;
			advance(it, 1);
			it->second = 100;
			advance(it, 1);
			it->second = 500;
			break;
		case 120:
			it->second = 300;
			advance(it, 1);
			it->second = 180;
			advance(it, 1);
			it->second = 200;
			advance(it, 1);
			it->second = 25;
			advance(it, 1);
			it->second = 60;
			advance(it, 1);
			it->second = 350;
			break;
		case 160:
			it->second = 300;
			advance(it, 1);
			it->second = 80;
			advance(it, 1);
			it->second = 100;
			advance(it, 1);
			it->second = 20;
			advance(it, 1);
			it->second = 20;
			advance(it, 1);
			it->second = 200;
			break;
		case 240:
			it->second = 150;
			advance(it, 1);
			it->second = 40;
			advance(it, 1);
			it->second = 50;
			advance(it, 1);
			it->second = 15;
			advance(it, 1);
			it->second = 20;
			advance(it, 1);
			it->second = UINT_MAX;
			break;
		}
	}
	else if (_caliber.first == "8U32")
	{
		it->second = 450;
		advance(it, 1);
		it->second = 300;
		advance(it, 1);
		it->second = 300;
		advance(it, 1);
		it->second = 60;
		advance(it, 1);
		it->second = 60;
		advance(it, 1);
		it->second = UINT_MAX;
	}
	else if (_caliber.first == "8U31")
	{
		it->second = 300;
		advance(it, 1);
		it->second = 100;
		advance(it, 1);
		it->second = 90;
		advance(it, 1);
		it->second = 25;
		advance(it, 1);
		it->second = 25;
		advance(it, 1);
		it->second = UINT_MAX;
	}
	else if (_caliber.first == "9K55")
	{
		it->second = 400;
		advance(it, 1);
		it->second = 150;
		advance(it, 1);
		it->second = 120;
		advance(it, 1);
		it->second = 25;
		advance(it, 1);
		it->second = 25;
		advance(it, 1);
		it->second = UINT_MAX;
	}
	else if (_caliber.first == "9K51")
	{
		it->second = 500;
		advance(it, 1);
		it->second = 240;
		advance(it, 1);
		it->second = 160;
		advance(it, 1);
		it->second = 35;
		advance(it, 1);
		it->second = 40;
		advance(it, 1);
		it->second = UINT_MAX;
	}
	else if (_caliber.first == "9K57")
	{
		if (_isCluster)
		{
			it->second = 40;
			advance(it, 1);
			it->second = 80;
			advance(it, 1);
			it->second = UINT_MAX;
			advance(it, 1);
			it->second = 1;
			advance(it, 1);
			it->second = 3;
		}
		else
		{
			it->second = 180;
			advance(it, 1);
			it->second = 120;
			advance(it, 1);
			it->second = 150;
			advance(it, 1);
			it->second = 7;
			advance(it, 1);
			it->second = 10;
		}

		advance(it, 1);
		it->second = UINT_MAX;
	}

	if (_caliber.first == "howitzer" || _caliber.first == "mortar")
	{
		if ((_range > 10.0))
		{
			for_each(targets.begin(), targets.end(), [&](auto& _n)
				{
					_n.second = static_cast<size_t> (_n.second * 1 / 10 * roundl(_range - 10.0));
				});
		}
		if (_isDoubleIndex) //чёт фигня какая-то с индексами ГРАУ, старые и новые, хз
		{
			it = targets.begin();
			advance(it, 1);
			it->second = static_cast<size_t>(it->second * 0.75);
			advance(it, 1);
			it->second = static_cast<size_t>(it->second * 0.75);
			advance(it, 1);
			it->second = static_cast<size_t>(it->second * 0.75);
			advance(it, 1);
			it->second = static_cast<size_t>(it->second * 0.75);
		}
		if (_isShortened)
		{
			for_each(targets.begin(), targets.end(), [&](auto& _n)
				{
					_n.second = static_cast<size_t> (_n.second * 1.5);
				});
		}
		
	}
	if (_isEntrenchedOrArmored)
	{
		it = targets.begin();
		advance(it, 1);
		it->second = static_cast<size_t>(it->second * 3); //РЛС, САУ, ещё чёто укрыты
		advance(it, 4);
		it->second = static_cast<size_t>(it->second * 3); //небронированная цель укрыта
	}
	else
	{
		it = targets.begin();
		it->second = static_cast<size_t>(it->second / 3); //артиллерия/минометы не укрыты
	}
	if (_isDestruction)
	{
		auto itTemp = targets.begin();
		advance(itTemp, 3);
		for_each(targets.begin(), itTemp, [&](auto& _n)
			{
				_n.second = static_cast<size_t> (_n.second * 3); //для уничтожения целей, указанных на подавление, число снарядов увеличивается втрое
			});
	}
	else
	{
		auto itTemp = targets.begin();
		advance(itTemp, 3);
		for_each(itTemp, targets.end(), [&](auto& _n)
			{
				_n.second = static_cast<size_t> (_n.second / 3); //всё наоборот, уменьшаем число снарядов втрое
			});
	}
	if (_isActiveReactive && (_caliber.second == 120 || _caliber.second == 152))
	{
		for_each(targets.begin(), targets.end(), [&](auto& _n)
			{
				_n.second = static_cast<size_t> (_n.second * 1.4); //всё наоборот, уменьшаем число снарядов втрое
			});
	}

	try
	{
		auto it2 = find_if(targets.begin(), targets.end(), [&](auto const& _n)
			{
				auto it3 = find_if(_n.first.begin(), _n.first.end(), [&](auto const& _m)
					{
						return _targetName == _m;
					});
				if (it3 != _n.first.end())
				{
					return true;
				}
				return false;

			});
		if (it2 == targets.end())
		{
			throw invalid_argument("Thiis target now found in the list for stationery targets");
		}
		else
		{
			return static_cast<size_t>(ceil(it2->second/ _numInBattery));
		}

	}
	catch (exception e)
	{
		cerr << e.what() << endl;
		exit(1);
	}
}
size_t Container::getNormMoving(pair<string, size_t> _caliber,bool _isCluster,  bool _isArmoured, size_t _numInBattery, string _targetName)
{
	map<vector<string>, size_t> targets = { {{"column"}, 0},
											{{"SPGArmoured", "SPMortarArmoured"}, 0 },
											{{"SPG", "SPMortar"}, 0},
											{{"launchSystem", "MLRS", "SAM", "AA", "helicopter"}, 0} };

	auto it = targets.begin();
	try
	{
		if (_caliber.first == "howitzer")
		{
			switch (_caliber.second)
			{
			case 76:
				it->second = 20;
				advance(it, 1);
				it->second = UINT_MAX;
				advance(it, 1);
				it->second = UINT_MAX;
				advance(it, 1);
				it->second = UINT_MAX;
				break;
			case 85:
				it->second = 16;
				advance(it, 1);
				it->second = UINT_MAX;
				advance(it, 1);
				it->second = UINT_MAX;
				advance(it, 1);
				it->second = UINT_MAX;
				break;
			case 100:
				it->second = 10;
				advance(it, 1);
				it->second = 16;
				advance(it, 1);
				it->second = 16;
				advance(it, 1);
				it->second = 18;
				break;
			case 120:
				it->second = 6;
				advance(it, 1);
				it->second = 12;
				advance(it, 1);
				it->second = 12;
				advance(it, 1);
				it->second = 9;
				break;
			case 122:
				it->second = 8;
				advance(it, 1);
				it->second = 16;
				advance(it, 1);
				it->second = 16;
				advance(it, 1);
				it->second = 10;
				break;
			case 130:
				it->second = 8;
				advance(it, 1);
				it->second = 10;
				advance(it, 1);
				it->second = 10;
				advance(it, 1);
				it->second = 8;
				break;
			case 152:

				it->second = 6;
				advance(it, 1);
				it->second = 10;
				advance(it, 1);
				it->second = 10;
				advance(it, 1);
				it->second = 8;
				break;
			case 203:
				throw runtime_error("Error! 203-mm howitzers can't be assigned to fire mission against moving targets!");
				break;
		}
	}
	else if (_caliber.first == "mortar")
	{
			throw runtime_error("Error! Mortars can't be assigned to fire mission against moving targets!");
	}
	else if (_caliber.first == "8U32")
	{
		it->second = 16;
		advance(it, 1);
		it->second = UINT_MAX;
		advance(it, 1);
		it->second = UINT_MAX;
		advance(it, 1);
		it->second = UINT_MAX;
	}
	else if (_caliber.first == "8U31")
	{
		it->second = 12;
		advance(it, 1);
		it->second = UINT_MAX;
		advance(it, 1);
		it->second = UINT_MAX;
		advance(it, 1);
		it->second = 12;
	}
	else if (_caliber.first == "9K55")
	{
		it->second = 36;
		advance(it, 1);
		it->second = 36;
		advance(it, 1);
		it->second = 36;
		advance(it, 1);
		it->second = 20;
	}
	else if (_caliber.first == "9K51")
		{
			it->second = 40;
			advance(it, 1);
			it->second = 40;
			advance(it, 1);
			it->second = 40;
			advance(it, 1);
			it->second = 30;
		}

	else if (_caliber.first == "9K57")
	{
		it->second = 16;
		advance(it, 1);
		it->second = UINT_MAX;
		advance(it, 1);
		if (_isCluster)
		{
			it->second = 10;			
		}
		else
		{
			it->second = 16;
		}
		advance(it, 1);
		it->second = 16;
		advance(it, 1);
		it->second = 3;
	}
	if (_isArmoured)
	{
		cout << "should be assigned at least one division!" << endl;
	}
	}
	catch (exception e)
	{
	cerr << e.what() << endl;
	exit(1);
	}
	try
	{
		auto it2 = find_if(targets.begin(), targets.end(), [&](auto const& _n)
			{
				auto it3 = find_if(_n.first.begin(), _n.first.end(), [&](auto const& _m)
					{
						return _targetName == _m;
					});
				if (it3 != _n.first.end())
				{
					return true;
				}
				return false;

			});
		if (it2 == targets.end())
		{
			throw invalid_argument("Thiis target now found in the list for stationery targets");
		}
		else
		{
			return static_cast<size_t>(ceil(it2->second/ _numInBattery));
		}

	}
	catch (exception e)
	{
		cerr << e.what() << endl;
		exit(1);
	}
}
size_t getNormFireAtWill(pair<string, size_t> _caliber, double _meters, double _minutes)
{
	size_t result = 0;
	try
	{
		if (_caliber.first == "howitzer")
		{
			switch (_caliber.second)
			{
			case 76: result = 18; break;
			case 85: result = 16; break;
			case 100: result = 12; break;
			case 120: result = 6; break;
			case 122: result = 8; break;
			case 130: result = 6; break;
			case 152: result = 6; break;
			case 203: throw runtime_error("this caliber is not used for fire-at-will missions"); break;
			}
		}
		else if (_caliber.first == "mortar")
		{
			switch (_caliber.second)
			{
			case 82: throw runtime_error("this caliber is not used for fire-at-will missions"); break;
			case 120: result = 6; break;
			case 160: result = 4; break;
			case 240: throw runtime_error("this caliber is not used for fire-at-will missions"); break;
			}
		}
		else
		{
			throw invalid_argument("MLRS cant be assigned to fire-at-will missions");
		}
	}
	catch (exception e)
	{
		cerr << e.what() << endl;
		exit(1);

	}

	return static_cast<size_t>(result * (_meters / 100) * (_minutes));
}
size_t getNormSmoke(pair<string, size_t> _caliber, double _meters, double _minutes, bool _isSideWind)
{
	size_t result = 0;
	try
	{
		if (_caliber.first == "howitzer" && _isSideWind)
		{
			switch (_caliber.second)
			{
			case 76: result = 4; break;
			case 85: result = 4; break;
			case 100: result = 2; break;
			case 120: result = 2; break;
			case 122: result = 1; break;
			case 130: throw runtime_error("this caliber is not used for fire-at-will missions"); break;
			case 152: result = 1; break;
			case 203: throw runtime_error("this caliber is not used for fire-at-will missions"); break;
			}
		}
		else if (_caliber.first == "howitzer" && not _isSideWind)
		{
			switch (_caliber.second)
			{
			case 76: result = 6; break;
			case 85: result = 6; break;
			case 100: result = 4; break;
			case 120: result = 2; break;
			case 122: result = 2; break;
			case 130: throw runtime_error("this caliber is not used for fire-at-will missions"); break;
			case 152: result = 2; break;
			case 203: throw runtime_error("this caliber is not used for fire-at-will missions"); break;
			}
		}
		else if (_caliber.first == "mortar" && _isSideWind)
		{
			switch (_caliber.second)
			{
			case 82: result = 4; break;
			case 120: result = 2; break;
			case 160: throw runtime_error("this caliber is not used for fire-at-will missions"); break;
			case 240: throw runtime_error("this caliber is not used for fire-at-will missions"); break;
			}
		}
		else if (_caliber.first == "mortar" && not _isSideWind)
		{
			switch (_caliber.second)
			{
			case 82: result = 6; break;
			case 120: result = 3; break;
			case 160: throw runtime_error("this caliber is not used for fire-at-will missions"); break;
			case 240: throw runtime_error("this caliber is not used for fire-at-will missions"); break;
			}
		}
		else
		{
			throw invalid_argument("MLRS cant be assigned to fire-at-will missions");
		}
	}
	catch (exception e)
	{
		cerr << e.what() << endl;
		exit(1);

	}
	return static_cast<size_t>(result * (_meters / 100) * (_minutes));
}
double Container::getRange(string _artyType)
{
	auto it = find_if(artyRangesAmmo.begin(), artyRangesAmmo.end(), [&](auto const& _n)
		{
			return _artyType == get<0>(_n);
		});

	return get<1>(*it);
}

double Container::getCaliber(string _artyType)
{
	auto it = find_if(artyRangesAmmo.begin(), artyRangesAmmo.end(), [&](auto const& _n)
		{
			return _artyType == get<0>(_n);
		});
	return get<4>(*it);
}

string Container::getType(string _artyType)
{
	try
	{
		list<string> howitzer = { "2C1", "2C3", "2C5", "2C7", "A-222", "2A31", "52P367", "2A18", "52P546", "52P482" };
		if (find(howitzer.begin(), howitzer.end(), _artyType) != howitzer.end()) { return "howitzer"; }

		list<string> mortar = { "2C4", "2C9", "2B9", "2B14", "52M853", "52M864", "2A51" };
		if (find(mortar.begin(), mortar.end(), _artyType) != mortar.end()) { return "mortar"; }

		list<string> MLRS = { "8U32", "8U32", "2K51", "2K55", "2K57" };
		if (find(MLRS.begin(), MLRS.end(), _artyType) != MLRS.end()) { return "MLRS"; }

		throw invalid_argument("Can't recognize such a type of arty!");
	}
	catch (exception e)
	{
		cerr << e.what() << endl;
		exit(1);

	}
}

bool Container::checkStationary(string _targetType)
{
	list<string> stat = { "artyEntranched", "mortarEntranched", "radarsAuto", "radioAuto", "AAOpened",
	"infEntranched", "CP", "tank", "APC", "infOpened", "autoCPOpened", "ATGMOpened", "ATGunOpened", "otherUnArmouredOpened" };
	if (find(stat.begin(), stat.end(), _targetType) != stat.end()) { return true; }
	return false;
}

size_t Container::measureFireRate(pair<string, size_t> _caliber, size_t _charge, double _temperature, size_t _numShells, string _specificType)
{
	map<size_t, size_t> targets = { {1, 0}, {3, 0 }, {5, 0}, {10, 0}, {15, 0}, {20, 0}, {25, 0}, {30, 0 }, {40, 0}, {50, 0},
											   {60, 0}, {120, 0}, {UINT_MAX, 0} };
	auto it = targets.begin();

	if (_caliber.first == "howitzer")
	{
		double percentUsageTemperature = 1.0;
		if (abs(20 - _temperature) > 10)
		{
			percentUsageTemperature += (_temperature - 20) / 100;
		}
		switch (_caliber.second)
		{
		case 76:
			it->second = static_cast<size_t>(ceil(15* percentUsageTemperature));
			advance(it, 1);
			it->second = static_cast<size_t>(ceil(35 * percentUsageTemperature));
			advance(it, 1);
			it->second = static_cast<size_t>(ceil(50 * percentUsageTemperature));
			advance(it, 1);
			it->second = static_cast<size_t>(ceil(70 * percentUsageTemperature));
			advance(it, 1);
			it->second = static_cast<size_t>(ceil(85 * percentUsageTemperature));
			advance(it, 1);
			it->second = static_cast<size_t>(ceil(100 * percentUsageTemperature));
			advance(it, 1);
			if (_charge == 0)
			{
	
			it->second = static_cast<size_t>(ceil(110 * percentUsageTemperature));
			advance(it, 1);
			it->second = static_cast<size_t>(ceil(115 * percentUsageTemperature));
			advance(it, 1);
			it->second = static_cast<size_t>(ceil(125 * percentUsageTemperature));
			advance(it, 1);
			it->second = static_cast<size_t>(ceil(140 * percentUsageTemperature));
			advance(it, 1);
			it->second = static_cast<size_t>(ceil(150 * percentUsageTemperature));
			advance(it, 1);
			it->second = static_cast<size_t>(ceil(220 * percentUsageTemperature));
			advance(it, 1);
			if (_numShells > 220) {
				size_t minutes = 120 + 60 * static_cast<size_t>((_numShells - 220) / 70) + static_cast<size_t>((_numShells - 220) % 70);
				
				it->second = static_cast<size_t>(ceil((220 + static_cast<size_t>(70 * (minutes - 120) / 60)) * percentUsageTemperature));
			}
			else { it->second = 220; }
			}
			else
			{
				
				it->second = static_cast<size_t>(ceil(115 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(130 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(160 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(180 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(200 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(320 * percentUsageTemperature));
				advance(it, 1);
				if (_numShells > 320) {
					size_t minutes = 120 + 60 * static_cast<size_t>((_numShells - 320) / 100) + static_cast<size_t>((_numShells - 320) % 70);

					it->second = static_cast<size_t>(ceil((320 + static_cast<size_t>(100 * (minutes - 120) / 60)) * percentUsageTemperature));
				}
				else { it->second = 220; }
				
			}
			break;
		case 85:
			it->second = static_cast<size_t>(ceil(10 * percentUsageTemperature));
			advance(it, 1);
			it->second = static_cast<size_t>(ceil(25 * percentUsageTemperature));
			advance(it, 1);
			it->second = static_cast<size_t>(ceil(40 * percentUsageTemperature));
			advance(it, 1);
			
			if (_charge == 0)
			{
				it->second = static_cast<size_t>(ceil(50 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(60 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(70 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(80 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(90 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(110 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(125 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(140 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(230 * percentUsageTemperature));
				advance(it, 1);
				if (_numShells > 230) {
					size_t minutes = 120 + 60 * static_cast<size_t>((_numShells - 230) / 80) + static_cast<size_t>((_numShells - 230) % 80);

					it->second = static_cast<size_t>(ceil((230 + static_cast<size_t>(80 * (minutes - 120) / 60)) * percentUsageTemperature));
				}
				else { it->second = 230; }
			}
			else
			{
				it->second = static_cast<size_t>(ceil(60 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(75 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(90 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(100 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(110 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(130 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(150 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(170 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(290 * percentUsageTemperature));
				advance(it, 1);
				if (_numShells > 290) {
					size_t minutes = 120 + 60 * static_cast<size_t>((_numShells - 290) / 100) + static_cast<size_t>((_numShells - 290) % 100);

					it->second = static_cast<size_t>(ceil((290 + static_cast<size_t>(100 * (minutes - 120) / 60)) * percentUsageTemperature));
				}
				else { it->second = 290; }
			}
			break;
		case 100:
			it->second = static_cast<size_t>(ceil(7 * percentUsageTemperature));
			advance(it, 1);
			it->second = static_cast<size_t>(ceil(18 * percentUsageTemperature));
			advance(it, 1);
			it->second = static_cast<size_t>(ceil(30 * percentUsageTemperature));
			advance(it, 1);
			it->second = static_cast<size_t>(ceil(50 * percentUsageTemperature));
			advance(it, 1);

			if (_charge == 0)
			{
				it->second = static_cast<size_t>(ceil(60 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(65 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(70 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(75 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(85 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(90 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(95 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(135 * percentUsageTemperature));
				advance(it, 1);
				if (_numShells > 135) {
					size_t minutes = 120 + 60 * static_cast<size_t>((_numShells - 135) / 40) + static_cast<size_t>((_numShells - 135) % 40);

					it->second = static_cast<size_t>(ceil((135 + static_cast<size_t>(40 * (minutes - 120) / 60)) * percentUsageTemperature));
				}
				else { it->second = 135; }
			}
			else
			{

				it->second = static_cast<size_t>(ceil(65 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(75 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(90 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(100 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(120 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(140 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(160 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(250 * percentUsageTemperature));
				advance(it, 1);
				if (_numShells > 250) {
					size_t minutes = 120 + 60 * static_cast<size_t>((_numShells - 250) / 80) + static_cast<size_t>((_numShells - 250) % 80);

					it->second = static_cast<size_t>(ceil((250 + static_cast<size_t>(80 * (minutes - 120) / 60)) * percentUsageTemperature));
				}
				else { it->second = 250; }
			}
			break;
		case 120:
			it->second = static_cast<size_t>(ceil(6 * percentUsageTemperature));
			advance(it, 1);
			it->second = static_cast<size_t>(ceil(16 * percentUsageTemperature));
			advance(it, 1);
			it->second = static_cast<size_t>(ceil(25 * percentUsageTemperature));
			advance(it, 1);
			it->second = static_cast<size_t>(ceil(40 * percentUsageTemperature));
			advance(it, 1);
			it->second = static_cast<size_t>(ceil(55 * percentUsageTemperature));
			advance(it, 1);
			if (_charge < 2)
			{
				it->second = static_cast<size_t>(ceil(65 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(70 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(75 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(85 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(90 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(100 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(150 * percentUsageTemperature));
				advance(it, 1);
				if (_numShells > 150) {
					size_t minutes = 120 + 60 * static_cast<size_t>((_numShells - 150) / 50) + static_cast<size_t>((_numShells - 150) % 50);

					it->second = static_cast<size_t>(ceil((150 + static_cast<size_t>(50 * (minutes - 120) / 60)) * percentUsageTemperature));
				}
				else { it->second = 150; }
			}
			else if (_charge > 3)
			{
				it->second = static_cast<size_t>(ceil(70 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(80 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(90 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(110 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(130 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(150 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(260 * percentUsageTemperature));
				advance(it, 1);
				if (_numShells > 260) {
					size_t minutes = 120 + 60 * static_cast<size_t>((_numShells - 260) / 80) + static_cast<size_t>((_numShells - 260) % 80);

					it->second = static_cast<size_t>(ceil((260 + static_cast<size_t>(80 * (minutes - 120) / 60)) * percentUsageTemperature));
				}
				else { it->second = 260; }
				
			}
			else
			{
				it->second = static_cast<size_t>(ceil(68 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(75 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(83 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(98 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(110 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(125 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(215 * percentUsageTemperature));
				advance(it, 1);
				if (_numShells > 215) {
					size_t minutes = 120 + 60 * static_cast<size_t>((_numShells - 215) / 65) + static_cast<size_t>((_numShells - 215) % 65);

					it->second = static_cast<size_t>(ceil((220 + static_cast<size_t>(65 * (minutes - 120) / 60)) * percentUsageTemperature));
				}
				else { it->second = 215; }

				}
			
			break;
		case 122:
			it->second = static_cast<size_t>(ceil(6 * percentUsageTemperature));
			advance(it, 1);
			it->second = static_cast<size_t>(ceil(16 * percentUsageTemperature));
			advance(it, 1);
			it->second = static_cast<size_t>(ceil(25 * percentUsageTemperature));
			advance(it, 1);
			it->second = static_cast<size_t>(ceil(40 * percentUsageTemperature));
			advance(it, 1);
			it->second = static_cast<size_t>(ceil(55 * percentUsageTemperature));
			advance(it, 1);
			if (_charge < 2)
			{
				it->second = static_cast<size_t>(ceil(65 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(70 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(75 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(85 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(90 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(100 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(150 * percentUsageTemperature));
				advance(it, 1);
				if (_numShells > 150) {
					size_t minutes = 120 + 60 * static_cast<size_t>((_numShells - 150) / 50) + static_cast<size_t>((_numShells - 150) % 50);

					it->second = static_cast<size_t>(ceil((150 + static_cast<size_t>(50 * (minutes - 120) / 60)) * percentUsageTemperature));
				}
				else { it->second = 150; }

				}
			
			else if (_charge > 3)
			{
				it->second = static_cast<size_t>(ceil(70 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(80 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(90 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(110 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(130 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(150 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(260 * percentUsageTemperature));
				advance(it, 1);
				if (_numShells > 260) {
					size_t minutes = 120 + 60 * static_cast<size_t>((_numShells - 260) / 80) + static_cast<size_t>((_numShells - 260) % 80);

					it->second = static_cast<size_t>(ceil((260 + static_cast<size_t>(80 * (minutes - 120) / 60)) * percentUsageTemperature));
				}
				else { it->second = 260; }
			}
			else
			{
				it->second = static_cast<size_t>(ceil(68 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(75 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(83 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(98 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(110 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(125 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(215 * percentUsageTemperature));
				advance(it, 1);
				if (_numShells > 215) {
					size_t minutes = 120 + 60 * static_cast<size_t>((_numShells - 215) / 65) + static_cast<size_t>((_numShells - 215) % 65);

					it->second = static_cast<size_t>(ceil((215 + static_cast<size_t>(65 * (minutes - 120) / 60)) * percentUsageTemperature));
				}
				else { it->second = 215; }
			}
			break;
		case 130:
			it->second = static_cast<size_t>(ceil(5 * percentUsageTemperature));
			advance(it, 1);
			it->second = static_cast<size_t>(ceil(12 * percentUsageTemperature));
			advance(it, 1);
			it->second = static_cast<size_t>(ceil(20 * percentUsageTemperature));
			advance(it, 1);
			it->second = static_cast<size_t>(ceil(35 * percentUsageTemperature));
			advance(it, 1);
			it->second = static_cast<size_t>(ceil(45 * percentUsageTemperature));
			advance(it, 1);
			it->second = static_cast<size_t>(ceil(55 * percentUsageTemperature));
			advance(it, 1);
			it->second = static_cast<size_t>(ceil(65 * percentUsageTemperature));
			advance(it, 1);

			if (_charge < 2)
			{
				it->second = static_cast<size_t>(ceil(70 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(80 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(90 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(100 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(160 * percentUsageTemperature));
				advance(it, 1);
				if (_numShells > 160) {
					size_t minutes = 120 + 60 * static_cast<size_t>((_numShells - 160) / 35) + static_cast<size_t>((_numShells - 160) % 35);

					it->second = static_cast<size_t>(ceil((160 + static_cast<size_t>(35 * (minutes - 120) / 60)) * percentUsageTemperature));
				}
				else { it->second = 160; }

			}
			else if (_charge > 2)
			{
				it->second = static_cast<size_t>(ceil(75 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(90 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(105 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(120 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(210 * percentUsageTemperature));
				advance(it, 1);
				if (_numShells > 210) {
					size_t minutes = 120 + 60 * static_cast<size_t>((_numShells - 210) / 70) + static_cast<size_t>((_numShells - 210) % 70);

					it->second = static_cast<size_t>(ceil((210 + static_cast<size_t>(70 * (minutes - 120) / 60)) * percentUsageTemperature));
				}
				else { it->second = 210; }

			}
			else
			{
				it->second = static_cast<size_t>(ceil(73 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(85 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(98 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(110 * percentUsageTemperature));
				advance(it, 1);
				it->second = static_cast<size_t>(ceil(185 * percentUsageTemperature));
				advance(it, 1);
				if (_numShells > 185) {
					size_t minutes = 120 + 60 * static_cast<size_t>((_numShells - 185) / 54) + static_cast<size_t>((_numShells - 185) % 54);

					it->second = static_cast<size_t>(ceil((185 + static_cast<size_t>(54 * (minutes - 120) / 60)) * percentUsageTemperature));
				}
				else { it->second = 185; }

			}
			break;
		case 152:
			it->second = static_cast<size_t>(ceil(4 * percentUsageTemperature));
			advance(it, 1);
			it->second = static_cast<size_t>(ceil(12 * percentUsageTemperature));
			advance(it, 1);
			it->second = static_cast<size_t>(ceil(20 * percentUsageTemperature));
			advance(it, 1);
			it->second = static_cast<size_t>(ceil(30 * percentUsageTemperature));
			advance(it, 1);
			it->second = static_cast<size_t>(ceil(40 * percentUsageTemperature));
			advance(it, 1);
			it->second = static_cast<size_t>(ceil(50 * percentUsageTemperature));
			advance(it, 1);
			if (_specificType == "howitzer")
			{
				it->second = static_cast<size_t>(ceil(60 * percentUsageTemperature));
				advance(it, 1);
				if (_charge < 2)
				{
					it->second = static_cast<size_t>(ceil(65 * percentUsageTemperature));
					advance(it, 1);
					it->second = static_cast<size_t>(ceil(75 * percentUsageTemperature));
					advance(it, 1);
					it->second = static_cast<size_t>(ceil(80 * percentUsageTemperature));
					advance(it, 1);
					it->second = static_cast<size_t>(ceil(90 * percentUsageTemperature));
					advance(it, 1);
					it->second = static_cast<size_t>(ceil(135 * percentUsageTemperature));
					advance(it, 1);
					if (_numShells > 135) {
						size_t minutes = 120 + 60 * static_cast<size_t>((_numShells - 135) / 45) + static_cast<size_t>((_numShells - 135) % 45);

						it->second = static_cast<size_t>(ceil((135 + static_cast<size_t>(45 * (minutes - 120) / 60)) * percentUsageTemperature));
					}
					else { it->second = 135; }

				}
				else if (_charge > 3)
				{
					it->second = static_cast<size_t>(ceil(70 * percentUsageTemperature));
					advance(it, 1);
					it->second = static_cast<size_t>(ceil(90 * percentUsageTemperature));
					advance(it, 1);
					it->second = static_cast<size_t>(ceil(105 * percentUsageTemperature));
					advance(it, 1);
					it->second = static_cast<size_t>(ceil(120 * percentUsageTemperature));
					advance(it, 1);
					it->second = static_cast<size_t>(ceil(210 * percentUsageTemperature));
					advance(it, 1);
					if (_numShells > 210) {
						size_t minutes = 120 + 60 * static_cast<size_t>((_numShells - 210) / 70) + static_cast<size_t>((_numShells - 210) % 70);

						it->second = static_cast<size_t>(ceil((210 + static_cast<size_t>(70 * (minutes - 120) / 60)) * percentUsageTemperature));
					}
					else { it->second = 210; }

				}
				else
				{
					it->second = static_cast<size_t>(ceil(68 * percentUsageTemperature));
					advance(it, 1);
					it->second = static_cast<size_t>(ceil(83 * percentUsageTemperature));
					advance(it, 1);
					it->second = static_cast<size_t>(ceil(93 * percentUsageTemperature));
					advance(it, 1);
					it->second = static_cast<size_t>(ceil(105 * percentUsageTemperature));
					advance(it, 1);
					it->second = static_cast<size_t>(ceil(173 * percentUsageTemperature));
					advance(it, 1);
					if (_numShells > 173) {
						size_t minutes = 120 + 60 * static_cast<size_t>((_numShells - 173) / 58) + static_cast<size_t>((_numShells - 173) % 58);

						it->second = static_cast<size_t>(ceil((173 + static_cast<size_t>(58 * (minutes - 120) / 60)) * percentUsageTemperature));
					}
					else { it->second = 173; }

				}

			}
			else if (_specificType == "howitzer-gun")
			{
				if (_charge < 2)
				{
					it->second = static_cast<size_t>(ceil(55 * percentUsageTemperature));
					advance(it, 1);
					it->second = static_cast<size_t>(ceil(60 * percentUsageTemperature));
					advance(it, 1);
					it->second = static_cast<size_t>(ceil(70 * percentUsageTemperature));
					advance(it, 1);
					it->second = static_cast<size_t>(ceil(75 * percentUsageTemperature));
					advance(it, 1);
					it->second = static_cast<size_t>(ceil(80 * percentUsageTemperature));
					advance(it, 1);
					it->second = static_cast<size_t>(ceil(120 * percentUsageTemperature));
					advance(it, 1);
					if (_numShells > 120) {
						size_t minutes = 120 + 60 * static_cast<size_t>((_numShells - 120) / 35) + static_cast<size_t>((_numShells - 120) % 35);

						it->second = static_cast<size_t>(ceil((220 + static_cast<size_t>(35 * (minutes - 120) / 60)) * percentUsageTemperature));
					}
					else { it->second = 120; }

				}
				else if (_charge > 5)
				{
					it->second = static_cast<size_t>(ceil(60 * percentUsageTemperature));
					advance(it, 1);
					it->second = static_cast<size_t>(ceil(70 * percentUsageTemperature));
					advance(it, 1);
					it->second = static_cast<size_t>(ceil(80 * percentUsageTemperature));
					advance(it, 1);
					it->second = static_cast<size_t>(ceil(95 * percentUsageTemperature));
					advance(it, 1);
					it->second = static_cast<size_t>(ceil(110 * percentUsageTemperature));
					advance(it, 1);
					it->second = static_cast<size_t>(ceil(200 * percentUsageTemperature));
					advance(it, 1);
					if (_numShells > 200) {
						size_t minutes = 120 + 60 * static_cast<size_t>((_numShells - 200) / 60) + static_cast<size_t>((_numShells - 200) % 60);

						it->second = static_cast<size_t>(ceil((200 + static_cast<size_t>(60 * (minutes - 120) / 60)) * percentUsageTemperature));
					}
					else { it->second = 200; }

				}
				else
				{
					it->second = static_cast<size_t>(ceil(58 * percentUsageTemperature));
					advance(it, 1);
					it->second = static_cast<size_t>(ceil(65 * percentUsageTemperature));
					advance(it, 1);
					it->second = static_cast<size_t>(ceil(75 * percentUsageTemperature));
					advance(it, 1);
					it->second = static_cast<size_t>(ceil(85 * percentUsageTemperature));
					advance(it, 1);
					it->second = static_cast<size_t>(ceil(95 * percentUsageTemperature));
					advance(it, 1);
					it->second = static_cast<size_t>(ceil(160 * percentUsageTemperature));
					advance(it, 1);
					if (_numShells > 160) {
						size_t minutes = 120 + 60 * static_cast<size_t>((_numShells - 160) / 48) + static_cast<size_t>((_numShells - 160) % 48);

						it->second = static_cast<size_t>(ceil((160 + static_cast<size_t>(48 * (minutes - 120) / 60)) * percentUsageTemperature));
					}
					else { it->second = 160; }
				}

			}
			else
			{
				it->second = static_cast<size_t>(ceil(60 * percentUsageTemperature));
				advance(it, 1);
				if (_charge < 2)
				{
					it->second = static_cast<size_t>(ceil(65 * percentUsageTemperature));
					advance(it, 1);
					it->second = static_cast<size_t>(ceil(75 * percentUsageTemperature));
					advance(it, 1);
					it->second = static_cast<size_t>(ceil(80 * percentUsageTemperature));
					advance(it, 1);
					it->second = static_cast<size_t>(ceil(90 * percentUsageTemperature));
					advance(it, 1);
					it->second = static_cast<size_t>(ceil(135 * percentUsageTemperature));
					advance(it, 1);
					if (_numShells > 135) {
						size_t minutes = 120 + 60 * static_cast<size_t>((_numShells - 135) / 45) + static_cast<size_t>((_numShells - 135) % 45);

						it->second = static_cast<size_t>(ceil((135 + static_cast<size_t>(45 * (minutes - 120) / 60)) * percentUsageTemperature));
					}
					else { it->second = 135; }
				}
				else if (_charge > 3)
				{
					it->second = static_cast<size_t>(ceil(70 * percentUsageTemperature));
					advance(it, 1);
					it->second = static_cast<size_t>(ceil(80 * percentUsageTemperature));
					advance(it, 1);
					it->second = static_cast<size_t>(ceil(95 * percentUsageTemperature));
					advance(it, 1);
					it->second = static_cast<size_t>(ceil(110 * percentUsageTemperature));
					advance(it, 1);
					it->second = static_cast<size_t>(ceil(200 * percentUsageTemperature));
					advance(it, 1);
					if (_numShells > 200) {
						size_t minutes = 120 + 60 * static_cast<size_t>((_numShells - 200) / 60) + static_cast<size_t>((_numShells - 200) % 60);

						it->second = static_cast<size_t>(ceil((200 + static_cast<size_t>(60 * (minutes - 120) / 60)) * percentUsageTemperature));
					}
					else { it->second = 200; }

				}
				else
				{
					it->second = static_cast<size_t>(ceil(68 * percentUsageTemperature));
					advance(it, 1);
					it->second = static_cast<size_t>(ceil(78 * percentUsageTemperature));
					advance(it, 1);
					it->second = static_cast<size_t>(ceil(88 * percentUsageTemperature));
					advance(it, 1);
					it->second = static_cast<size_t>(ceil(100 * percentUsageTemperature));
					advance(it, 1);
					it->second = static_cast<size_t>(ceil(168 * percentUsageTemperature));
					advance(it, 1);
					if (_numShells > 168) {
						size_t minutes = 120 + 60 * static_cast<size_t>((_numShells - 168) / 53) + static_cast<size_t>((_numShells - 168) % 53);

						it->second = static_cast<size_t>(ceil((168 + static_cast<size_t>(53 * (minutes - 120) / 60)) * percentUsageTemperature));
					}
					else { it->second = 168; }

				}
			}
			break;
		case 203:
			it->second = static_cast<size_t>(ceil(1 * percentUsageTemperature));
			advance(it, 1);
			it->second = static_cast<size_t>(ceil(2 * percentUsageTemperature));
			advance(it, 1);
			it->second = static_cast<size_t>(ceil(3 * percentUsageTemperature));
			advance(it, 1);
			it->second = static_cast<size_t>(ceil(6 * percentUsageTemperature));
			advance(it, 1);
			it->second = static_cast<size_t>(ceil(9 * percentUsageTemperature));
			advance(it, 1);
			it->second = static_cast<size_t>(ceil(12 * percentUsageTemperature));
			advance(it, 1);
			it->second = static_cast<size_t>(ceil(15 * percentUsageTemperature));
			advance(it, 1);
			it->second = static_cast<size_t>(ceil(18 * percentUsageTemperature));
			advance(it, 1);
			it->second = static_cast<size_t>(ceil(22 * percentUsageTemperature));
			advance(it, 1);
			it->second = static_cast<size_t>(ceil(26 * percentUsageTemperature));
			advance(it, 1);
			it->second = static_cast<size_t>(ceil(30 * percentUsageTemperature));
			advance(it, 1);
			it->second = static_cast<size_t>(ceil(60 * percentUsageTemperature));
			advance(it, 1);
			if (_numShells > 60) {
				size_t minutes = 120 + 60 * static_cast<size_t>((_numShells - 60) / 25) + static_cast<size_t>((_numShells - 60) % 25);

				it->second = static_cast<size_t>(ceil((60 + static_cast<size_t>(25 * (minutes - 120) / 60)) * percentUsageTemperature));
			}
			else { it->second = 60; }
			break;
		}
	}
	else if (_caliber.first == "mortar")
	{
		switch (_caliber.second)
		{
		case 82:
			it->second = 20;
			advance(it, 1);
			it->second = 45;
			advance(it, 1);
			
			advance(it, 1);
			if (_charge == 0)
			{
				it->second = 75;
				advance(it, 1);
				it->second = 110;
				advance(it, 1);
				it->second = 125;
				it->second = 65;
				advance(it, 1);
				it->second = 150;
				advance(it, 1);
				it->second = UINT_MAX;
				advance(it, 1);
				if (_numShells > 150) {

					it->second = UINT_MAX;
				}
				else { it->second = 150; }
			}
			else if (_charge == 1)
			{
				it->second = 60;
				advance(it, 1);
				it->second = 75;
				advance(it, 1);
				it->second = 85;
				it->second = 100;
				advance(it, 1);
				it->second = 140;
				advance(it, 1);
				it->second = UINT_MAX;
				advance(it, 1);
				if (_numShells > 140) {
					
					it->second = UINT_MAX;
				}
				else { it->second = 140; }
				
			}
			break;
		case 120:
			it->second = 9;
			advance(it, 1);
			it->second = 25;
			advance(it, 1); 
			it->second = 30;
			advance(it, 1);
			it->second = 35;
			advance(it, 1);
			it->second = 40;
			advance(it, 1);

			if (_charge == 0)
			{
				it->second = 50;
				advance(it, 1);
				it->second = 70;
				advance(it, 1);
				it->second = 110;
				advance(it, 1);
				if (_numShells > 110) {
					it->second = 150;
				}
				else
				{
					it->second = 110;
				}
			}
			else if (_charge > 0)
			{
				it->second = 75;
				advance(it, 1);
				it->second = 105;
				advance(it, 1);
				it->second = 165;
				advance(it, 1);
				if (_numShells > 165) {
					it->second = 225;
				}
				else
				{
					it->second = 165;
				}
			}
			break;
		case 160:
			it->second = 3;
			advance(it, 1);
			it->second = 9;
			advance(it, 1);
			it->second = 12;
			advance(it, 1);
			it->second = 18;
			advance(it, 1);
			it->second = 21;
			advance(it, 1);
			it->second = 30;
			advance(it, 1);
			it->second = 48;
			advance(it, 1);
			it->second = 88;
			if (_numShells > 88) {
				it->second = 128;
			}
			else
			{
				it->second = 88;
			}
			break;
		case 240:
			it->second = 1;
			advance(it, 1);
			it->second = 3;
			advance(it, 1);
			it->second = 5;
			advance(it, 1);
			it->second = 10;
			advance(it, 1);

			if (_charge == 0)
			{
				it->second = 15;
				advance(it, 1);
				it->second = 30;
				advance(it, 1);
				it->second = 38;
				advance(it, 1);
				it->second = 51;
				advance(it, 1);
				if (_numShells > 51) {
					it->second = 64;
				}
				else
				{
					it->second = 51;
				}
			}
			else if (_charge > 0)
			{
				it->second = 13;
				advance(it, 1);
				it->second = 20;
				advance(it, 1);
				it->second = 25;
				advance(it, 1);
				it->second = 35;
				advance(it, 1);

				if (_numShells > 35) {
					it->second = 45;
				}
				else
				{
					it->second = 35;
				}

			}
			break;
		}
	}
	else if (_caliber.first == "8U32") //здесь вообще по=хорошему надо считать отдельно нормы на батарею и на дивизион
	{
			it->second = UINT_MAX;
			advance(it, 1);
			it->second = 16;
			advance(it, 1);
			it->second = 16;
			advance(it, 1);
			it->second = 16*2;
			advance(it, 1);
			it->second = 16*3;
			advance(it, 1);
			it->second = 16*3;
			advance(it, 1);
			it->second = 16 * 4;
			advance(it, 1);
			it->second = 16 * 5;
			advance(it, 1);
			it->second = 16 * 6;
			advance(it, 1);
			it->second = 16 * 7;
			advance(it, 1);
			it->second = 16 * 8;
			advance(it, 1);
			if (_numShells > 16*8) {
				size_t minutes = 120 + 60 * static_cast<size_t>((_numShells - 16*8) / 16*6) + static_cast<size_t>((_numShells - 16*8) % 16*6);

				it->second = 16*8 + static_cast<size_t>(16*6 * (minutes - 120) / 60);
			}
			else { it->second = 16*8; }
	}
	else if (_caliber.first == "8U31")
	{
		it->second = UINT_MAX; //1
		advance(it, 1);
		it->second = 12; //3
		advance(it, 1);
		it->second = 12; //5
		advance(it, 1);
		it->second = 12; //10
		advance(it, 1);
		it->second = 12 * 2; //15
		advance(it, 1);
		it->second = 12 * 2; //20
		advance(it, 1);
		it->second = 12 * 2; //25
		advance(it, 1);
		it->second = 12 * 3; //30
		advance(it, 1);
		it->second = 12 * 3; //40
		advance(it, 1);
		it->second = 12 * 4; //50
		advance(it, 1);
		it->second = 12 * 5; //60
		advance(it, 1);
		if (_numShells > 12*5) {
			size_t minutes = 120 + 60 * static_cast<size_t>((_numShells - 12*5) / 12*3) + static_cast<size_t>((_numShells - 12*5) % 12*3);

			it->second = 12*5 + static_cast<size_t>(12*3 * (minutes - 120) / 60);
		}
		else { it->second = 12*5; }
	}
	else if (_caliber.first == "9K55")
	{
		it->second = UINT_MAX; //1
		advance(it, 1);
		it->second = 36; //3
		advance(it, 1);
		it->second = 36 * 2; //5
		advance(it, 1);
		it->second = 36 * 2; //10
		advance(it, 1);
		it->second = 36 * 3; //15
		advance(it, 1);
		it->second = 36 * 3; //20
		advance(it, 1);
		it->second = 36 * 4; //25
		advance(it, 1);
		it->second = 36 * 5; //30
		advance(it, 1);
		it->second = 36 * 6; //40
		advance(it, 1);
		it->second = 36 * 7; //50
		advance(it, 1);
		it->second = 36 * 8; //60
		advance(it, 1);
		if (_numShells > 36*8) {
			size_t minutes = 120 + 60 * static_cast<size_t>((_numShells - 36*8) / 36*6) + static_cast<size_t>((_numShells - 36*8) % 36*6);

			it->second = 36*8 + static_cast<size_t>(36*6 * (minutes - 120) / 60);
		}
		else { it->second = 36*8; }
	}
	else if (_caliber.first == "9K51")
	{
		it->second = UINT_MAX; //1
		advance(it, 1);
		it->second = 40; //3
		advance(it, 1);
		it->second = 40 * 2; //5
		advance(it, 1);
		it->second = 40 * 2; //10
		advance(it, 1);
		it->second = 40 * 2; //15
		advance(it, 1);
		it->second = 40 * 3; //20
		advance(it, 1);
		it->second = 40 * 3; //25
		advance(it, 1);
		it->second = 40 * 4; //30
		advance(it, 1);
		it->second = 40 * 4; //40
		advance(it, 1);
		it->second = 40 * 5; //50
		advance(it, 1);
		it->second = 40 * 6; //60
		advance(it, 1);
		if (_numShells > 40*6) {
			size_t minutes = 120 + 60 * static_cast<size_t>((_numShells - 40*6) / 40*4) + static_cast<size_t>((_numShells - 40*6) % 40*4);

			it->second = 40*6 + static_cast<size_t>(40*4 * (minutes - 120) / 60);
		}
		else { it->second = 40*6; }
	}
	else if (_caliber.first == "9K57")
	{
		it->second = UINT_MAX; //1
		advance(it, 1);
		it->second = 40; //3
		advance(it, 1);
		it->second = 40 * 2; //5
		advance(it, 1);
		it->second = 40 * 2; //10
		advance(it, 1);
		it->second = 40 * 2; //15
		advance(it, 1);
		it->second = 40 * 3; //20
		advance(it, 1);
		it->second = 40 * 3; //25
		advance(it, 1);
		it->second = 40 * 4; //30
		advance(it, 1);
		it->second = 40 * 4; //40
		advance(it, 1);
		it->second = 40 * 5; //50
		advance(it, 1);
		it->second = 40 * 6; //60
		advance(it, 1);
		if (_numShells > 40*6) {
			size_t minutes = 120 + 60 * static_cast<size_t>((_numShells - 40 * 6) / 40 * 4) + static_cast<size_t>((_numShells - 40*6) % 40 * 4);

			it->second = 40*6 + static_cast<size_t>(40 * 4 * (minutes - 120) / 60);
		}
		else { it->second = 40 * 6; }

	}
	try
	{
		auto it2 = find_if(targets.begin(), targets.end(), [&](auto const& _n)
			{
				return _n.second > _numShells; //возвращаем норму по потолку
			});
		auto it3 = it2;
		if (it2 != targets.begin())
		{
			advance(it3, -1);
		}
		

		if (it2->second - it3->second == 0 && it2 != targets.begin())
		{
			throw runtime_error("Could't interpolate: equal neighbouring shell numbers");
		}
		else if (it2->second - it3->second == 0 && it2 == targets.begin())
		{
			return it3->first;
		}

		if (it2 != targets.end())
		{
			double percent = (_numShells - it3->second) / (it2->second - it3->second);
			return it3->first + static_cast<size_t>(percent * (it2->first - it3->first)); //и интерполируем эту величину 
		}
		else 
		{
			throw invalid_argument("Can't find fitting thresholds for this number of shells");
		}

	}
	catch (exception e)
	{
		cerr << e.what() << endl;
		exit(1);
	}
}

string Container::getSpecificType(string _artyType)
{
	try
	{
		if (getType(_artyType) != "howitzer")
		{
			throw invalid_argument("Trying to implement specification of howitzers on mortars or MLRS");
		}
	}
	catch (invalid_argument e)
	{
		cerr << e.what() << endl;
		exit(1);
	}

	//такое себе, как оказалось, это было даже излишним. Между ними разницы немного

	list<string> howitzer = { "2C1","2A31","2A18"};
	if (find(howitzer.begin(), howitzer.end(), _artyType) != howitzer.end()) { return "howitzer"; }

	list<string> gunhowitzer = {"2C3", "2C5", "2C7", "A-222","52P367","52P546"};
	if (find(gunhowitzer.begin(), gunhowitzer.end(), _artyType) != gunhowitzer.end()) { return "howitzer-gun"; }

	list<string> howitzergun = {"52P482" };
	if (find(howitzergun.begin(), howitzergun.end(), _artyType) != howitzergun.end()) { return "gun-howitzer"; }
}

size_t getCharge(string _artyName, string _targetName, string _fuzeType, bool _isEntrenchedOrArmored, double _range)
{
	if (_artyName == "2A18")
	{
		if (_isEntrenchedOrArmored) return 0;

		if (_range < 5.800) return 4;
		else if (_range < 7.000) return 3;
		else if (_range < 8.200) return 2;
		else if (_range < 9.400) return 1;
		else return 0;
	}

	if (_artyName == "52P546")
	{
		if (_isEntrenchedOrArmored) return 0;

		else if (_range < 6.000) return 6;
		else if (_range < 7.000) return 5;
		else if (_range < 9.000) return 4;
		else if (_range < 10.200) return 3;
		else if (_range < 11.400) return 2;
		else if (_range < 13.600) return 1;
		else return 0;
	}

	return 0;
}