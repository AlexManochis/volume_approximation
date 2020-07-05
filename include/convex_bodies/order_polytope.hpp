// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.
//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef ORDER_POLYTOPE
#define ORDER_POLYTOPE

#include <limits>
#include <iostream>
#include "solve_lp.h"

//min and max values for the Hit and Run functions


// H-polytope class
template <typename Point>
class OrderPolytope{
public:
    typedef Point PointType;
    typedef typename Point::FT NT;
    typedef typename std::vector<NT>::iterator viterator;
    //using RowMatrixXd = Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    //typedef RowMatrixXd MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;

private:
    MT E; //matrix E
    MT A;
    VT b; // vector b, s.t.: Ax<=b
    unsigned int            _d; //dimension
    std::pair<Point,NT> _inner_ball;
    //NT maxNT = 1.79769e+308;
    //NT minNT = -1.79769e+308;
    NT sqd;
    NT maxNT = std::numeric_limits<NT>::max();
    NT minNT = std::numeric_limits<NT>::lowest();

public:

    OrderPolytope()
    {
        typedef typename Point::FT NT;
        _inner_ball = ComputeChebychevBall<NT, Point>(A, b);
    }

    std::pair<Point,NT> InnerBall() const
    {
        return _inner_ball;
    }

    //Compute Chebyshev ball of H-polytope P:= Ax<=b
    //Use LpSolve library
    std::pair<Point,NT> ComputeInnerBall()
    {
        _inner_ball = ComputeChebychevBall<NT, Point>(A, b);
        return _inner_ball;
    }


    // return dimension
    unsigned int dimension() const {
        return _d;
    }


    // return the number of facets
    int num_of_hyperplanes() const {
        return A.rows();
    }

    int num_of_generators() const {
        return 0;
    }


    // return the matrix A
    MT get_mat() const {
        return A;
    }


    // return the vector b
    VT get_vec() const {
        return b;
    }

    MT get_AA() const {
        return A * A.transpose();
    }

    // change the matrix A
    void set_mat(const MT &A2) {
        A = A2;
    }


    // change the vector b
    void set_vec(const VT &b2) {
        b = b2;
    }

    Point get_mean_of_vertices() const {
        return Point(_d);
    }

    NT get_max_vert_norm() const {
        return 0.0;
    }

    void init(const unsigned int dim, const MT &_E, const MT &_A) {
        _d = dim;
        E = _E;
        A = _A;
        b = VT::Zero(2 * dim + _E.rows());
        for (unsigned int i = 0; i < 2 * dim; ++i) {
            b(i) = 1.0;
        }

        NT row_norm;
        for (int i = 0; i < num_of_hyperplanes(); ++i) {
            row_norm = A.row(i).norm();
            A.row(i) = A.row(i) / row_norm;
            b(i) = b(i) / row_norm;
        }

        //std::cout<<A<<"\n"<<std::endl;
        //std::cout<<E<<"\n"<<std::endl;
        //std::cout<<b<<"\n"<<std::endl;
        sqd = std::sqrt(2.0);
    }

    //define matrix A and vector b, s.t. Ax<=b and the dimension
    void init(const std::vector<std::vector<NT> > &Pin) {
        _d = Pin[0][1] - 1;
        //A.resize(Pin.size() - 1, _d);
        b.resize(Pin.size() - 1);
        for (unsigned int i = 1; i < Pin.size(); i++) {
            b(i - 1) = Pin[i][0];
            for (unsigned int j = 1; j < _d + 1; j++) {
                //A(i - 1, j - 1) = -Pin[i][j];
            }
        }
    }


    // print polytope in input format
    void print() {
        std::cout << " " << E.rows() << " " << _d << " float" << std::endl;
        for (unsigned int i = 0; i < E.rows(); i++) {
            for (unsigned int j = 0; j < _d; j++) {
                std::cout << E(i, j) << " ";
            }
            std::cout << "<= " << b(i) << std::endl;
        }
    }


    // Compute the reduced row echelon form
    // used to transofm {Ax=b,x>=0} to {A'x'<=b'}
    // e.g. Birkhoff polytopes
    /*
    // Todo: change the implementation in order to use eigen matrix and vector.
    int rref(){
        to_reduced_row_echelon_form(_A);
        std::vector<int> zeros(_d+1,0);
        std::vector<int> ones(_d+1,0);
        std::vector<int> zerorow(_A.size(),0);
        for (int i = 0; i < _A.size(); ++i)
        {
            for (int j = 0; j < _d+1; ++j){
                if ( _A[i][j] == double(0)){
                    ++zeros[j];
                    ++zerorow[i];
                }
                if ( _A[i][j] == double(1)){
                    ++ones[j];
                }
            }
        }
        for(typename stdMatrix::iterator mit=_A.begin(); mit<_A.end(); ++mit){
            int j =0;
            for(typename stdCoeffs::iterator lit=mit->begin(); lit<mit->end() ; ){
                if(zeros[j]==_A.size()-1 && ones[j]==1)
                    (*mit).erase(lit);
                else{ //reverse sign in all but the first column
                    if(lit!=mit->end()-1) *lit = (-1)*(*lit);
                    ++lit;
                }
                ++j;
            }
        }
        //swap last and first columns
        for(typename stdMatrix::iterator mit=_A.begin(); mit<_A.end(); ++mit){
            double temp=*(mit->begin());
            *(mit->begin())=*(mit->end()-1);
            *(mit->end()-1)=temp;
        }
        //delete zero rows
        for (typename stdMatrix::iterator mit=_A.begin(); mit<_A.end(); ){
            int zero=0;
            for(typename stdCoeffs::iterator lit=mit->begin(); lit<mit->end() ; ++lit){
                if(*lit==double(0)) ++zero;
            }
            if(zero==(*mit).size())
                _A.erase(mit);
            else
                ++mit;
        }
        //update _d
        _d=(_A[0]).size();
        // add unit vectors
        for(int i=1;i<_d;++i){
            std::vector<double> e(_d,0);
            e[i]=1;
            _A.push_back(e);
        }
        // _d should equals the dimension
        _d=_d-1;
        return 1;
    }*/


    //Check if Point p is in H-polytope P:= Ax<=b
    int is_in(const Point &p) const
    {
        const NT* b_data = b.data();

        for (int i = 0; i < _d; i++) {
            //Check if corresponding hyperplane is violated
            if (*b_data - p[i] < NT(0))
                return 0;

            b_data++;
        }
        for (int i = 0; i < _d; i++) {
            //Check if corresponding hyperplane is violated
            if (*b_data + p[i] < NT(0))
                return 0;

            b_data++;
        }
        for (int i = 0; i < E.rows(); i++) {
            //Check if corresponding hyperplane is violated
            if (*b_data + (p[E(i,0)-1] - p[E(i,1)-1])/sqd < NT(0))
                return 0;

            b_data++;
        }
        return -1;
    }


    // compute intersection point of ray starting from r and pointing to v
    // with polytope discribed by A and b
    std::pair<NT,NT> line_intersect(Point &r, Point &v) const
    {

        NT lamda = 0, min_plus = NT(maxNT), max_minus = NT(minNT);
        VT sum_nom, sum_denom;
        //unsigned int i, j;
        unsigned int j;
        int m = num_of_hyperplanes();


        sum_nom.noalias() = b - A * r.getCoefficients();
        sum_denom.noalias() = A * v.getCoefficients();

        NT* sum_nom_data = sum_nom.data();
        NT* sum_denom_data = sum_denom.data();

        for (int i = 0; i < m; i++) {

            if (*sum_denom_data == NT(0)) {
                //std::cout<<"div0"<<std::endl;
                ;
            } else {
                lamda = *sum_nom_data / *sum_denom_data;
                if (lamda < min_plus && lamda > 0) min_plus = lamda;
                if (lamda > max_minus && lamda < 0) max_minus = lamda;
            }

            sum_nom_data++;
            sum_denom_data++;
        }
        return std::pair<NT, NT>(min_plus, max_minus);
    }


    // compute intersection points of a ray starting from r and pointing to v
    // with polytope discribed by A and b
    std::pair<NT,NT> line_intersect(Point &r, Point &v, VT& Ar,
                                    VT& Av, bool pos = false) const {


        NT lamda = 0, min_plus = NT(maxNT), max_minus = NT(minNT);
        VT sum_nom;
        int m = num_of_hyperplanes(), facet;

        Ar.noalias() = A * r.getCoefficients();
        sum_nom = b - Ar;
        Av.noalias() = A * v.getCoefficients();;


        NT* Av_data = Av.data();
        NT* sum_nom_data = sum_nom.data();

        for (int i = 0; i < m; i++) {
            if (*Av_data == NT(0)) {
                //std::cout<<"div0"<<std::endl;
                ;
            } else {
                lamda = *sum_nom_data / *Av_data;
                if (lamda < min_plus && lamda > 0) {
                    min_plus = lamda;
                    if (pos) facet = i;
                }else if (lamda > max_minus && lamda < 0) max_minus = lamda;
            }

            Av_data++;
            sum_nom_data++;
        }
        if (pos) return std::pair<NT, NT>(min_plus, facet);
        return std::pair<NT, NT>(min_plus, max_minus);
    }

    std::pair<NT,NT> line_intersect(Point &r, Point &v, VT& Ar,
                                    VT& Av, const NT &lambda_prev, bool pos = false) const
    {
        NT lamda = 0, min_plus = NT(maxNT), max_minus = NT(minNT);
        VT sum_nom;
        NT mult;
        //unsigned int i, j;
        unsigned int j;
        int m = num_of_hyperplanes(), facet;

        Ar.noalias() += lambda_prev*Av;
        sum_nom = b - Ar;
        Av.noalias() = A * v.getCoefficients();

        NT* sum_nom_data = sum_nom.data();
        NT* Av_data = Av.data();

        for (int i = 0; i < m; i++) {
            if (*Av_data == NT(0)) {
                //std::cout<<"div0"<<std::endl;
                ;
            } else {
                lamda = *sum_nom_data / *Av_data;
                if (lamda < min_plus && lamda > 0) {
                    min_plus = lamda;
                    if (pos) facet = i;
                }else if (lamda > max_minus && lamda < 0) max_minus = lamda;
            }
            Av_data++;
            sum_nom_data++;
        }
        if (pos) return std::pair<NT, NT>(min_plus, facet);
        return std::pair<NT, NT>(min_plus, max_minus);
    }


    // compute intersection point of a ray starting from r and pointing to v
    // with polytope discribed by A and b
    std::pair<NT, int> line_positive_intersect(Point &r,
                                               Point &v,
                                               VT& Ar,
                                               VT& Av) const {
        return line_intersect(r, v, Ar, Av, true);
    }


    // compute intersection point of a ray starting from r and pointing to v
    // with polytope discribed by A and b
    std::pair<NT, int> line_positive_intersect(Point &r,
                                               Point &v,
                                               VT& Ar,
                                               VT& Av,
                                               const NT &lambda_prev) const
    {
        return line_intersect(r, v, Ar, Av, lambda_prev, true);
    }

    //------------------------------------------------------------------------------------
    // compute intersection point of a ray starting from r and pointing to v
    // with polytope discribed by A and b
    template <typename update_parameters>
    std::pair<NT, int> line_first_positive_intersect(Point const& r,
                                                     Point const& v,
                                                     VT& Ar,
                                                     VT& Av,
                                                     update_parameters &params) const
    {
        NT lamda = 0, min_plus = NT(maxNT);
        VT sum_nom;
        //unsigned int i, j;
        unsigned int j;
        int m = num_of_hyperplanes(), facet;
        //viterator rit, vit, Ariter = Ar.begin(), Aviter = Av.begin();

        NT* Ar_data = Ar.data();
        NT* Av_data = Av.data();

        for (int i = 0; i < _d; ++i) {
            *Ar_data = r[i];
            *Av_data = v[i];
            Av_data++;
            Ar_data++;
        }
        for (int i = 0; i < _d; ++i) {
            *Ar_data = -r[i];
            *Av_data = -v[i];
            Av_data++;
            Ar_data++;
        }
        for (int i = 0; i < E.rows(); ++i) {
            *Ar_data = (-r[E(i,0)-1] + r[E(i,1)-1])/sqd;
            *Av_data = (-v[E(i,0)-1] + v[E(i,1)-1])/sqd;
            Av_data++;
            Ar_data++;
        }
        //std::cout<<v.getCoefficients().transpose()<<"\n"<<std::endl;
        //std::cout<<(Ar.transpose()).transpose()<<"\n"<<std::endl;
        //std::cout<<A*r.getCoefficients()<<"\n"<<std::endl;
        //std::cout<<Av.transpose()<<"\n"<<std::endl;
        //std::cout<<(A*v.getCoefficients()).transpose()<<"\n"<<std::endl;

        //Ar.noalias() = A * r.getCoefficients();
        sum_nom.noalias() = b - Ar;
        //Av.noalias() = A * v.getCoefficients();

        Av_data = Av.data();
        NT* sum_nom_data = sum_nom.data();

        for (int i = 0; i < m; i++) {
            if (*Av_data == NT(0)) {
                //std::cout<<"div0"<<std::endl;
                ;
            } else {
                lamda = *sum_nom_data / *Av_data;
                if (lamda < min_plus && lamda > 0) {
                    min_plus = lamda;
                    //if (pos){
                    facet = i;
                    //facet_k = i;
                    params.inner_vi_ak = *Av_data;
                    //}
                }
            }

            Av_data++;
            sum_nom_data++;
        }
        params.facet_prev = facet;
        return std::pair<NT, int>(min_plus, facet);
        //return line_intersect(r, v, Ar, Av, lambda_prev, true);
    }


    template <typename update_parameters>
    std::pair<NT, int> line_positive_intersect(Point const& r,
                                               Point const& v,
                                               VT& Ar,
                                               VT& Av,
                                               NT const& lambda_prev,
                                               MT const& AA,
                                               update_parameters &params) const
    {
        NT lamda = 0, min_plus = NT(maxNT);
        VT sum_nom;//
        // , sum_denom, sum2;
        //unsigned int i, j;
        unsigned int j;
        int m = num_of_hyperplanes(), facet;
        NT inner_prev = params.inner_vi_ak;
        //viterator vit, Ariter = Ar.begin(), Aviter = Av.begin();

        //std::cout<<"[3]facet_k = "<<facet_k<<", inner_vi_ak = "<<inner_prev<<std::endl;
        Ar.noalias() += lambda_prev*Av;
        if(params.hit_ball) {
            Av.noalias() += (-2.0 * inner_prev) * (Ar / params.ball_inner_norm);
        } else {
            Av.noalias() += (-2.0 * inner_prev) * AA.col(params.facet_prev);
            //sum2 = (-2.0 * inner_prev) * ((*Ariter)/params.ball_inner_norm);
        }
        sum_nom.noalias() = b - Ar;

        NT* sum_nom_data = sum_nom.data();
        NT* Av_data = Av.data();

        for (int i = 0; i < m; i++) {
            if (*Av_data == NT(0)) {
                //std::cout<<"div0"<<std::endl;
                ;
            } else {
                lamda = *sum_nom_data / *Av_data;
                if (lamda < min_plus && lamda > 0) {
                    min_plus = lamda;
                    facet = i;
                    //facet_k = i;
                    params.inner_vi_ak = *Av_data;
                }
            }
            Av_data++;
            sum_nom_data++;
        }
        params.facet_prev = facet;
        return std::pair<NT, int>(min_plus, facet);
    }


    template <typename update_parameters>
    std::pair<NT, int> line_positive_intersect(Point const& r,
                                               Point const& v,
                                               VT& Ar,
                                               VT& Av,
                                               NT const& lambda_prev,
                                               update_parameters &params) const
    {
        NT lamda = 0, min_plus = NT(maxNT);
        VT sum_nom;
        unsigned int j;
        int m = num_of_hyperplanes(), facet;

        Ar.noalias() += lambda_prev * Av;
        sum_nom.noalias() = b - Ar;

        NT* Av_data = Av.data();

        for (int i = 0; i < _d; ++i) {
            *Av_data = v[i];
            Av_data++;
        }
        for (int i = 0; i < _d; ++i) {
            *Av_data = -v[i];
            Av_data++;
        }
        for (int i = 0; i < E.rows(); ++i) {
            *Av_data = (-v[E(i, 0)-1] + v[E(i, 1)-1])/sqd;
            Av_data++;
        }

        //Av.noalias() = A * v.getCoefficients();

        NT* sum_nom_data = sum_nom.data();
        Av_data = Av.data();

        for (int i = 0; i < m; i++) {
            if (*Av_data == NT(0)) {
                //std::cout<<"div0"<<std::endl;
                ;
            } else {
                lamda = *sum_nom_data / *Av_data;
                if (lamda < min_plus && lamda > 0) {
                    min_plus = lamda;
                    facet = i;
                    //facet_k = i;
                    params.inner_vi_ak = *Av_data;
                }
            }
            Av_data++;
            sum_nom_data++;
        }
        params.facet_prev = facet;
        return std::pair<NT, int>(min_plus, facet);
    }

    //-----------------------------------------------------------------------------------//

    //First coordinate ray intersecting convex polytope
    std::pair<NT,NT> line_intersect_coord(Point &r,
                                          const unsigned int &rand_coord,
                                          VT& lamdas) const
    {

        NT lamda = 0, min_plus = NT(maxNT), max_minus = NT(minNT);
        VT sum_denom;
        int m = num_of_hyperplanes();

        sum_denom = A.col(rand_coord);
        lamdas.noalias() = b - A * r.getCoefficients();

        NT* lamda_data = lamdas.data();
        NT* sum_denom_data = sum_denom.data();

        for (int i = 0; i < m; i++) {

            if (*sum_denom_data == NT(0)) {
                //std::cout<<"div0"<<sum_denom<<std::endl;
                ;
            } else {
                lamda = *lamda_data * (1 / *sum_denom_data);
                if (lamda < min_plus && lamda > 0) min_plus = lamda;
                if (lamda > max_minus && lamda < 0) max_minus = lamda;

            }
            lamda_data++;
            sum_denom_data++;
        }
        return std::pair<NT, NT>(min_plus, max_minus);
    }


    //Not the first coordinate ray intersecting convex
    std::pair<NT,NT> line_intersect_coord(Point &r,
                                          const Point &r_prev,
                                          const unsigned int rand_coord,
                                          const unsigned int rand_coord_prev,
                                          VT& lamdas) const
    {
        NT lamda = 0, min_plus = NT(maxNT), max_minus = NT(minNT);

        int m = num_of_hyperplanes();

        lamdas.noalias() += A.col(rand_coord_prev)* (r_prev[rand_coord_prev] - r[rand_coord_prev]);
        NT* data = lamdas.data();

        for (int i = 0; i < m; i++) {
            NT a = A(i, rand_coord);

            if (a == NT(0)) {
                //std::cout<<"div0"<<std::endl;
                ;
            } else {
                lamda = *data / a;
                if (lamda < min_plus && lamda > 0) min_plus = lamda;
                if (lamda > max_minus && lamda < 0) max_minus = lamda;

            }
            data++;
        }
        return std::pair<NT, NT>(min_plus, max_minus);
    }


    // Apply linear transformation, of square matrix T^{-1}, in H-polytope P:= Ax<=b
    void linear_transformIt(const MT &T) {
        A = A * T;
    }


    // shift polytope by a point c
    void shift(const VT &c)
    {
        b -= A*c;
    }


    // return for each facet the distance from the origin
    std::vector<NT> get_dists(const NT &radius) const
    {
        unsigned int i=0;
        std::vector <NT> dists(num_of_hyperplanes(), NT(0));
        typename std::vector<NT>::iterator disit = dists.begin();
        for ( ; disit!=dists.end(); disit++, i++)
            *disit = b(i) / A.row(i).norm();

        return dists;
    }


    // no points given for the rounding, you have to sample from the polytope
    template <typename T>
    bool get_points_for_rounding (const T &randPoints) {
        return false;
    }

    MT get_T() const {
        return A;
    }

    void normalize() {

        //NT row_norm;
        //for (int i = 0; i < num_of_hyperplanes(); ++i) {
        //    row_norm = A.row(i).norm();
        //    A.row(i) = A.row(i) / row_norm;
        //    b(i) = b(i) / row_norm;
        //}

    }

    void compute_reflection(Point &v, const Point &, const int facet) const
    {
        v += -2 * v.dot(A.row(facet)) * A.row(facet);
    }

    //void compute_reflection(Point &v, const Point &, const NT &inner_vi_ak, const int &facet) const {

    //Point a((-2.0 * inner_vi_ak) * A.row(facet));
    //v += a;
    //}

    template <typename update_parameters>
    void compute_reflection(Point &v, const Point &, update_parameters const& params) const {

        Point s = v;
        //Point q(A.row(params.facet_prev));
        if (params.facet_prev < _d) {
            //std::cout << "hi1 " << std::endl;
            //s.set_coord(params.facet_prev, -2.0 * params.inner_vi_ak);
            v.set_coord(params.facet_prev, -2.0 * params.inner_vi_ak + v[params.facet_prev]);
        } else if (params.facet_prev < 2 * _d) {
            //std::cout << "hi2 " << std::endl;
           // std::cout<<"inner_vec = "<<params.inner_vi_ak<<"\n"<<std::endl;
            //std::cout<<"v = "<<v.getCoefficients().transpose()<<"\n"<<std::endl;
            //s.set_coord(params.facet_prev - _d, 2.0 * params.inner_vi_ak);
            //std::cout<<"A.row = "<<A.row(params.facet_prev)<<std::endl;
            v.set_coord(params.facet_prev - _d, 2.0 * params.inner_vi_ak + v[params.facet_prev - _d]);
        } else {
           // std::cout << "hi3 " << std::endl;
            v.set_coord(E(params.facet_prev - 2*_d, 0) - 1, (-2.0 * params.inner_vi_ak) *
                            A(params.facet_prev, E(params.facet_prev - 2*_d, 0) - 1) + v[E(params.facet_prev - 2*_d, 0) - 1]);
            v.set_coord(E(params.facet_prev - 2*_d, 1) - 1, (-2.0 * params.inner_vi_ak) *
                            A(params.facet_prev, E(params.facet_prev - 2*_d, 1) - 1) + v[E(params.facet_prev - 2*_d, 1) - 1]);

            //s.set_coord(E(params.facet_prev - 2 * _d, 0) - 1, (-2.0 * params.inner_vi_ak) * s[E(params.facet_prev - 2 * _d, 0) - 1]);
            //s.set_coord(E(params.facet_prev - 2 * _d, 1) - 1, (-2.0 * params.inner_vi_ak) * s[E(params.facet_prev - 2 * _d, 1) - 1]);
            //s *= (NT(1)/(sqd*sqd));
        }
        Point a((-2.0 * params.inner_vi_ak) * A.row(params.facet_prev));
        //std::cout<<"a = "<<a.getCoefficients().transpose()<<"\n"<<std::endl;

        s += a;
        //v=s;
        //std::cout<<"s = "<<s.getCoefficients().transpose()<<"\n"<<std::endl;
        //std::cout<<"v = "<<v.getCoefficients().transpose()<<"\n"<<std::endl;
        if (!(v==s)) {
            std::cout << "v != s " << std::endl;
        }
    }

    void free_them_all() {}

};

#endif
