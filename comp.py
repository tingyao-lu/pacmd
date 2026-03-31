from openmm import NoseHooverChain


# thermostat
# NoseHooverChain(double, double, double, double, int, int, int, int, int, vector<int>, vector<pair<int,int>>)
x = NoseHooverChain.get(300, 10, 10, 10, 1, 3, 5, 2, 3, [9], [(9, 1)])
