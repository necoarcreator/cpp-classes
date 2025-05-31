void Test1(size_t N_x, vector<double>* p, vector<double> *rho, vector<double> *v){
    double p1, rho1, v1;
    double p2, rho2, v2;
    rho1 = 1.0;
    v1 = 0.0;
    p1 = 1.0;
    rho2 = 0.125;
    v2 = 0.0;
    p2 = 0.1;
    for(int i = 0; i < N_x; ++i){
        if (i < 0.5 * N_x){
            p[i] = p1;
            rho[i] = rho1;
            v[i] = v1;
        }else{
            p[i] = p2;
            rho[i] = rho2;
            v[i] = v2;
        }
    }
}
void Test2(size_t N_x, vector<double>* p, vector<double> *rho, vector<double> *v){
    double p1, rho1, v1;
    double p2, rho2, v2;
    rho1 = 1.0;
    v1 = -2.0;
    p1 = 0.4;
    rho2 = 1.0;
    v2 = 2.0;
    p2 = 0.4;
    for(int i = 0; i < N_x; ++i){
        if (i < 0.5 * N_x){
            p[i] = p1;
            rho[i] = rho1;
            v[i] = v1;
        }else{
            p[i] = p2;
            rho[i] = rho2;
            v[i] = v2;
        }
    }
}
void Test3(size_t N_x, vector<double>* p, vector<double> *rho, vector<double> *v){
    double p1, rho1, v1;
    double p2, rho2, v2;
    rho1 = 1.0;
    v1 = 0.0;
    p1 = 1000.0;
    rho2 = 1.0;
    v2 = 0.0;
    p2 = 0.01;
    for(int i = 0; i < N_x; ++i){
        if (i < 0.5 * N_x){
            p[i] = p1;
            rho[i] = rho1;
            v[i] = v1;
        }else{
            p[i] = p2;
            rho[i] = rho2;
            v[i] = v2;
        }
    }
}
void Test4(size_t N_x, vector<double>* p, vector<double> *rho, vector<double> *v){
    double p1, rho1, v1;
    double p2, rho2, v2;
    rho1 = 1.0;
    v1 = 0.0;
    p1 = 0.01;
    rho2 = 1.0;
    v2 = 0.0;
    p2 = 100.0;
    for(int i = 0; i < N_x; ++i){
        if (i < 0.5 * N_x){
            p[i] = p1;
            rho[i] = rho1;
            v[i] = v1;
        }else{
            p[i] = p2;
            rho[i] = rho2;
            v[i] = v2;
        }
    }
}
void Test5(size_t N_x, vector<double>* p, vector<double> *rho, vector<double> *v){
    double p1, rho1, v1;
    double p2, rho2, v2;
    rho1 = 5.99924;
    v1 = 19.5975;
    p1 = 460.894;
    rho2 = 5.99242;
    v2 = -6.19633;
    p2 = 46.0950;
    for(int i = 0; i < N_x; ++i){
        if (i < 0.5 * N_x){
            p[i] = p1;
            rho[i] = rho1;
            v[i] = v1;
        }else{
            p[i] = p2;
            rho[i] = rho2;
            v[i] = v2;
        }
    }
}

