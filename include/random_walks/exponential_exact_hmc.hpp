// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

// Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_WALKS_EXPONENTIAL_EXACT_HMC_WALK_HPP
#define RANDOM_WALKS_EXPONENTIAL_EXACT_HMC_WALK_HPP

#include "convex_bodies/ball.h"
#include "convex_bodies/ballintersectconvex.h"
#include "convex_bodies/hpolytope.h"
#include "convex_bodies/vpolytope.h"
#include "convex_bodies/vpolyintersectvpoly.h"
#include "convex_bodies/zpolytope.h"
#include "convex_bodies/zonoIntersecthpoly.h"


// Billiard walk which accelarates each step for uniform distribution

struct hmc_exponential
{
    hmc_exponential(double L)
            :   param(L, true)
    {}

    hmc_exponential()
            :   param(0, false)
    {}

    struct parameters
    {
        parameters(double L, bool set)
                :   m_L(L), set_L(set)
        {}
        double m_L;
        bool set_L;
    };

    struct update_parameters
    {
        update_parameters()
                :   facet_prev(0), inner_vi_ak(0.0)
        {}
        int facet_prev;
        double inner_vi_ak;
    };

    parameters param;


    template
            <
                    typename Polytope,
                    typename RandomNumberGenerator
            >
    struct Walk
    {
        typedef typename Polytope::PointType Point;
        typedef typename Polytope::MT MT;
        typedef typename Polytope::VT VT;
        typedef typename Point::FT NT;

        template <typename GenericPolytope>
        Walk(GenericPolytope const& P, Point const& p, VT const& c, NT const& Temp, RandomNumberGenerator &rng)
        {
            _update_params = update_parameters();
           _L = compute_diameter<GenericPolytope>
                ::template compute<NT>(P);
            _AA.noalias()= P.get_AA();
            _c = c;
            _Ac = P.get_mat() * _c;
            _Temp = Temp;
            initialize(P, p, rng);
        }

        template <typename GenericPolytope>
        Walk(GenericPolytope const& P, Point const& p, VT const& c, NT const& Temp, RandomNumberGenerator &rng,
             parameters const& params)
        {
            _update_params = update_parameters();
            _L = params.set_L ? params.m_L
                              : compute_diameter<GenericPolytope>
                                ::template compute<NT>(P);
            _AA.noalias()= P.get_AA();
            _c = c;
            _Ac = P.get_mat() * _c;
            _Temp = Temp;
            initialize(P, p, rng);
        }

        template
                <
                        typename GenericPolytope
                >
        inline void apply(GenericPolytope const& P,
                          Point &p,   // a point to start
                          unsigned int const& walk_length,
                          RandomNumberGenerator &rng)
        {
            unsigned int n = P.dimension();
            NT T = rng.sample_urdist() * _L;
            const NT dl = 0.995;
            int it;

            for (auto j=0u; j<walk_length; ++j)
            {
                T = rng.sample_urdist() * _L;
                _v = GetDirection<Point>::apply(n, rng);
                Point p0 = _p;

                it = 0;
                std::pair<NT, int> pbpair = P.quadratic_positive_intersect(_p, _v, _Ac, _c, _Temp, _lambdas, _Av,
                                                                           _lambda_prev, _update_params);
                std::cout<<"it = "<<it<<", t = "<<pbpair.first<<", T ="<<T<<std::endl;
                if (T <= pbpair.first) {
                    _p += ((T * _v) - (((T*T)/(2*_Temp)) * _c));
                    std::cout<<"_v = "<<_v.getCoefficients().transpose()<<"\n"<<std::endl;
                    std::cout<<"_p = "<<_p.getCoefficients().transpose()<<"\n"<<std::endl;
                    _lambda_prev = T;
                    continue;
                }

                _lambda_prev = dl * pbpair.first;
                _p += ((_lambda_prev * _v) - (((_lambda_prev*_lambda_prev)/(2.0*_Temp)) * _c));
                T -= _lambda_prev;
                P.compute_reflection(_v, _c, _Temp, _Ac, _lambda_prev,_update_params);
                it++;

                while (it < 100*n)
                {
                    std::pair<NT, int> pbpair
                            = P.quadratic_positive_intersect(_p, _v, _Ac, _c, _Temp, _lambdas, _Av, _lambda_prev, _AA,
                                                             _update_params);
                    std::cout<<"it = "<<it<<", t = "<<pbpair.first<<", T ="<<T<<std::endl;
                    if (T <= pbpair.first) {
                        _p += ((T * _v) - (((T*T)/(2.0*_Temp)) * _c));
                        std::cout<<"_v = "<<_v.getCoefficients().transpose()<<"\n"<<std::endl;
                        std::cout<<"_p = "<<_p.getCoefficients().transpose()<<"\n"<<std::endl;
                        _lambda_prev = T;
                        break;
                    }
                    _lambda_prev = dl * pbpair.first;
                    _p += ((_lambda_prev * _v) - (((_lambda_prev*_lambda_prev)/(2.0*_Temp)) * _c));
                    T -= _lambda_prev;
                    P.compute_reflection(_v, _c, _Temp, _Ac, _lambda_prev, _update_params);
                    it++;
                }
                if (it == 100*n){
                    std::cout<<"limit reached"<<std::endl;
                    _p = p0;
                }
            }
            p = _p;
        }

        inline void update_delta(NT L)
        {
            _L = L;
        }

    private :

        template
                <
                        typename GenericPolytope
                >
        inline void initialize(GenericPolytope const& P,
                               Point const& p,
                               RandomNumberGenerator &rng)
        {
            //std::cout<<
            unsigned int n = P.dimension();
            const NT dl = 0.995;
            _lambdas.setZero(P.num_of_hyperplanes());
            _Av.setZero(P.num_of_hyperplanes());
            _p = p;
            _v = GetDirection<Point>::apply(n, rng);

            NT T = rng.sample_urdist() * _L;
            Point p0 = _p;
            int it = 0;
            std::cout<<"_c = "<<_c.transpose()<<", Temperature = "<<_Temp<<std::endl;

            std::pair<NT, int> pbpair
                    = P.quadratic_first_positive_intersect(_p, _v, _Ac, _c, _Temp, _lambdas, _Av, _update_params);
            std::cout<<"t first = "<<pbpair.first<<", T ="<<T<<std::endl;
            std::cout<<"is_in = "<<P.is_in(_p + (((pbpair.first*0.9) * _v) - ((pbpair.first*0.9*pbpair.first*0.9)/(2.0*_Temp)) * _c));
            if (T <= pbpair.first) {
                _p += ((T * _v) - ((T*T)/(2.0*_Temp)) * _c);
                std::cout<<"_v = "<<_v.getCoefficients().transpose()<<"\n"<<std::endl;
                std::cout<<"_p = "<<_p.getCoefficients().transpose()<<"\n"<<std::endl;
                _lambda_prev = T;
                return;
            }
            _lambda_prev = dl * pbpair.first;
            _p += ((_lambda_prev * _v) - ((_lambda_prev*_lambda_prev)/(2.0*_Temp)) * _c);
            T -= _lambda_prev;
            P.compute_reflection(_v, _c, _Temp, _Ac, _lambda_prev, _update_params);

            while (it < 100*n)
            {
                std::pair<NT, int> pbpair
                        = P.quadratic_positive_intersect(_p, _v, _Ac, _c, _Temp, _lambdas, _Av, _lambda_prev, 
                                                         _AA, _update_params);
                std::cout<<"it = "<<it<<", t = "<<pbpair.first<<", T ="<<T<<std::endl;
                if (T <= pbpair.first) {
                    _p += ((T * _v) - ((T*T)/(2.0*_Temp)) * _c);
                    std::cout<<"vec_increase = "<<((T * _v) - ((T*T)/(2.0*_Temp)) * _c).getCoefficients().transpose()<<"\n"<<std::endl;
                    std::cout<<"_p = "<<_p.getCoefficients().transpose()<<"\n"<<std::endl;
                    _lambda_prev = T;
                    break;
                } else if (it == 100*n) {
                    _lambda_prev = rng.sample_urdist() * pbpair.first;
                    _p += ((_lambda_prev * _v) - ((_lambda_prev*_lambda_prev)/(2.0*_Temp)) * _c);
                    break;
                }
                _lambda_prev = dl * pbpair.first;
                _p += ((_lambda_prev * _v) - ((_lambda_prev*_lambda_prev)/(2.0*_Temp)) * _c);
                T -= _lambda_prev;
                P.compute_reflection(_v, _c, _Temp, _Ac, _lambda_prev, _update_params);
                it++;
            }
        }

        double _L;
        Point _p;
        Point _v;
        NT _lambda_prev, _Temp;
        MT _AA;
        update_parameters _update_params;
        VT _lambdas;
        VT _Av;
        VT _c;
        VT _Ac;
    };

};


#endif


