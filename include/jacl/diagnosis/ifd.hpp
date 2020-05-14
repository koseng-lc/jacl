/**
*   @author : koseng (Lintang)
*   @brief : Instrument Fault Detection
*/

#pragma once

#include <boost/variant.hpp>

#include <jacl/traits.hpp>
#include <jacl/system/continuous.hpp>
#include <jacl/system/discrete.hpp>
#include <jacl/pole_placement.hpp>
#include <jacl/lti_common.hpp>
#include <jacl/plotter.hpp>
#include <jacl/medium.hpp>

namespace jacl{ namespace diagnosis{

template <typename _System>
class IFD:public ::jacl::system::detail::BaseSystemClient<typename _System::base_t>{
public:
    using diag_pack_t = std::array<std::tuple<arma::vec, arma::vec, bool>, _System::n_outputs>;

public:    
    IFD(_System* _sys, std::initializer_list<double> _threshold)
        : sys_(_sys)
        , aux_sys_(nullptr)
        , plt_(&aux_sys_, plotSig<_System::n_outputs-1>::idx(), sys_->dt()){  
        
        auto _th_it = _threshold.begin();
        for(int i(0); i < threshold_.size() && _th_it != _threshold.end(); i++){
            threshold_[i] = *_th_it++;
        }
     
        initDOs<_System::n_outputs-1, decltype(dos_)>::init(&dos_, sys_);
        setSubjectDOs<_System::n_outputs-1, decltype(dos_)>::set(&dos_, sys_);
    }
    ~IFD(){}
    virtual void init(std::initializer_list<
                        typename arma::Col<
                            typename _System::scalar_t>::template fixed<_System::n_states-1>> _poles,
                      std::string&& _plot_title = "",
                      std::initializer_list<std::string> _plot_name = {}){
        std::vector<typename arma::Col<
                        typename _System::scalar_t>::template fixed<_System::n_states-1>> poles(_poles);
        setPoleDOs<_System::n_outputs-1, decltype(dos_)>::set(&dos_, poles, ifd_tag());
        if(_plot_name.size()){
            plt_.init();
            plt_.setTitle(std::move(_plot_title));
            plt_.setPlotName(_plot_name);
        }        
    }
    virtual auto detect(const typename _System::input_t& _in,
                        const typename _System::output_t& _out)
                        -> diag_pack_t{
        diag_pack_t res;
        std::vector<typename _System::output_t> y_hat;
        typename _System::output_t delta_y;
        getEst<_System::n_outputs-1, decltype(dos_)>::get(&dos_, _in, _out, &y_hat);
        for(int i(0); i < _System::n_outputs; i++){
            delta_y = _out - y_hat[i];
            auto err(1.);
            for(int j(0); j < _System::n_outputs; j++){
                if(j != i)
                    err *= delta_y(j);
            }
            std::get<0>(res[i]) = y_hat[i];
            std::get<1>(res[i]) = delta_y;
            std::get<2>(res[i]) = err > threshold_[i];
        }
        //-- temporary yet - it's not the actual implementation
        //-- we're gonna have many estimation signals
        aux_sys_.collectSig(y_hat[0], _in, y_hat[0]);
        return res;
    }
    auto viewSignals(){
        if(plt_.isInitialized())
            plt_.start();
    }
protected:
    template <typename __System = _System, std::size_t chosen_state=0>
    class DedicatedObserver:public ::jacl::system::BaseSystem<typename __System::state_space_t>{
    protected:
        using A11_t = typename arma::Mat<typename __System::scalar_t>::template fixed<1,1>;
        using A12_t = typename arma::Mat<typename __System::scalar_t>::template fixed<1,_System::n_states-1>;
        using A21_t = typename arma::Mat<typename __System::scalar_t>::template fixed<_System::n_states-1,1>;
        using A22_t = typename arma::Mat<typename __System::scalar_t>::template fixed<_System::n_states-1,_System::n_states-1>;
        using B1_t = typename arma::Mat<typename __System::scalar_t>::template fixed<1,_System::n_inputs>;
        using B2_t = typename arma::Mat<typename __System::scalar_t>::template fixed<_System::n_states-1,_System::n_inputs>;

    public:
        static constexpr double CHOSEN_STATE{chosen_state};
        DedicatedObserver(__System* _sys=nullptr)
            : ::jacl::system::BaseSystem<typename __System::state_space_t>(
                    ::jacl::system::detail::BaseSystemClient<typename __System::base_t>::ss(_sys))
            , z_( arma::zeros<arma::vec>(__System::n_outputs-1, 1) )
            , prev_z_(z_){
        }
        ~DedicatedObserver(){}
         //-- overload
        auto convolve(const typename __System::input_t& _in, const typename __System::output_t& _meas){
            meas_ = _meas;
            typename __System::output_t est( convolve(_in) );
            return est;
        }    
        auto setPole(const typename arma::Col<
                            typename __System::scalar_t>::template fixed<__System::n_states-1>& _poles){
            poles_ = _poles;
            updateVar();
        }
    protected:
        auto convolve(const typename __System::input_t& _in)
            -> typename __System::output_t override{
            setIn(_in);
            this->state_ = dstate();
            //-- it must be not here, just adapt with the reality
            this->prev_state_ = this->state_;
            //--
            this->out_ = output();            
            return this->out_;
        }
        void setIn(const typename __System::input_t& _in){
            this->in_ = _in;
        }
        template <typename T = __System>
        auto dstate(const arma::mat& term1,
                    const arma::mat& term2,
                    const arma::mat& term3,
                    const arma::vec& _yI)
            -> typename std::enable_if_t<
                ::jacl::traits::is_continuous_system<T>::value, typename __System::state_t>{
            //-- not implemented yet
            return arma::vec(T::n_outputs, 1, arma::fill::zeros);
        }
        template <typename T = __System>
        auto dstate(const arma::mat& _term1,
                    const arma::mat& _term2,
                    const arma::mat& _term3,
                    const arma::vec& _yI)
            -> typename std::enable_if_t<
                ::jacl::traits::is_discrete_system<T>::value,
                typename arma::Col<typename __System::scalar_t>::template fixed<__System::n_states-1>>{
            return _term1*prev_z_ + _term1*L_*_yI
                    + _term2*_yI + _term3*this->in_; 
        }
        auto dstate() -> typename __System::state_t override{
            arma::mat term1( A22_ - L_*A12_ );
            arma::mat term2( A21_ - L_*A11_ );
            arma::mat term3( B2_ - L_*B1_ );
            arma::vec yI( meas_(CHOSEN_STATE) - this->ss_->D().row(CHOSEN_STATE) * this->in_ );
            z_ = dstate(term1, term2, term3, yI);
            q_hat_ = arma::join_vert(yI, prev_z_ + L_*yI);
            prev_z_ = z_;
            return  arma::inv(M_) * q_hat_;
        }
        auto output() -> typename __System::output_t override{
            return this->ss_->C() * this->prev_state_ + this->ss_->D() * this->in_;
        }
        void updateVar() override{
            transform();
            if(!An_.is_empty() && !Bn_.is_empty() && !poles_.is_empty()){
                A11_ = A11(An_);
                A12_ = A12(An_);
                A21_ = A21(An_);
                A22_ = A22(An_);
                B1_ = B1(Bn_);
                B2_ = B2(Bn_);
                jacl::state_space::Linear<double, _System::n_states-1,1,1> ss(A22_, arma::zeros(_System::n_states-1,1),
                                                                    A12_, arma::zeros(1,1));
                jacl::pole_placement::KautskyNichols(&ss, poles_, &L_, jacl::pole_placement::PolePlacementType::Observer);
                L_.print("DOS GAIN : ");
                // if(jacl::lti_common::observable(ss.A(), ss.C()))
                //     std::cout << "SS Observable !" << std::endl;
                // arma::cx_mat eigvec;
                // arma::cx_vec eigval;
                // arma::eig_gen(eigval, eigvec, ss.A() - gain*ss.C());
                // eigval.print("EigVal : ");
            }
        }
        auto transform(){
            if(!this->ss_->A().is_empty()
                && !this->ss_->B().is_empty()
                && !this->ss_->C().is_empty()){
                arma::mat C_t = arma::trans( this->ss_->C() );
                arma::mat M_t(__System::n_states, __System::n_states, arma::fill::eye);
                for(int j(0); j < __System::n_states; j++){
                    if(C_t(j, CHOSEN_STATE) != 0){
                        M_t = arma::shift(M_t, -j, 1);
                        break;
                    }
                }
                for(int j(0); j < __System::n_states; j++){
                    if(C_t(j,CHOSEN_STATE) != 0){
                        M_t(j,0) = C_t(j,CHOSEN_STATE);
                    }
                }
                M_ = arma::trans(M_t);
                An_ = M_ * this->ss_->A() * arma::inv(M_);
                Bn_ = M_ * this->ss_->B(); 
            }                         
        }
        inline auto A11(const typename __System::state_space_t::state_matrix_t& _A)
            -> A11_t{
            return _A.submat(0,0,0,0);
        }
        inline auto A12(const typename __System::state_space_t::state_matrix_t& _A)
            -> A12_t{
            return _A.submat(0,1,0,__System::n_states-1);
        }
        inline auto A21(const typename __System::state_space_t::state_matrix_t& _A)
            -> A21_t{
            return _A.submat(1,0,__System::n_states-1,0);
        }
        inline auto A22(const typename __System::state_space_t::state_matrix_t& _A)
            -> A22_t{
            return _A.submat(1,1,__System::n_states-1,__System::n_states-1);
        }
        inline auto B1(const typename __System::state_space_t::input_matrix_t& _B)
            -> B1_t{
            return _B.head_rows(1);
        }
        inline auto B2(const typename __System::state_space_t::input_matrix_t& _B)
            -> B2_t{
            return _B.tail_rows(__System::n_states-1);
        }
    protected:
        typename __System::state_space_t::state_matrix_t An_;
        typename __System::state_space_t::input_matrix_t Bn_;

        A11_t A11_;
        A12_t A12_;
        A21_t A21_;
        A22_t A22_;
        B1_t B1_;
        B2_t B2_;
        typename __System::state_space_t::state_matrix_t M_;
        typename arma::Col<typename __System::scalar_t>::template fixed<__System::n_states-1> L_;
        arma::vec z_;
        arma::vec prev_z_;
        typename __System::state_t q_hat_;
        typename __System::output_t meas_;
        arma::vec poles_;
    };

protected:
    std::array<double, _System::n_outputs> threshold_;
    _System* sys_;

    template <int n, typename T>
    struct DOVariant;
    template <int n, typename ...Args>
    struct DOVariant<n, boost::variant<Args...>>{
        using type = typename DOVariant<n-1, boost::variant<Args..., DedicatedObserver<_System,n-1>>>::type;
    };
    template <typename ...Args>
    struct DOVariant<1, boost::variant<Args...>>{
        using type = typename boost::variant<Args..., DedicatedObserver<_System,0>>;
    };
    std::vector<typename DOVariant<_System::n_outputs-1,
                boost::variant<DedicatedObserver<_System, _System::n_outputs-1>>>::type> dos_;

    template <typename _StateSpace>
    class AuxSys:public ::jacl::system::BaseSystem<_StateSpace>{
    public:
        AuxSys(_StateSpace* _ss)
            : ::jacl::system::BaseSystem<_StateSpace>(_ss){}
        auto collectSig(const typename ::jacl::system::BaseSystem<_StateSpace>::state_t& _state,
                        const typename ::jacl::system::BaseSystem<_StateSpace>::input_t& _in,
                        const typename ::jacl::system::BaseSystem<_StateSpace>::output_t& _out){
            this->prev_state_ = _state;
            this->in_ = _in;
            this->out_ = _out;            
        }        
    protected:
        void setIn(const typename ::jacl::system::BaseSystem<_StateSpace>::input_t& _in){}
        auto dstate() -> typename ::jacl::system::BaseSystem<_StateSpace>::state_t{}
        auto output() -> typename ::jacl::system::BaseSystem<_StateSpace>::output_t{}
    };
    //-- so we gonna have same dimension in all signals
    AuxSys<typename _System::state_space_t> aux_sys_;
    jacl::Plotter<AuxSys<typename _System::state_space_t>> plt_;

    template <int n, typename DOs>
    struct initDOs{
        static constexpr std::size_t CHOSEN_STATE{(_System::n_outputs-1) - n};
        static auto init(DOs* _dos, _System* _sys){
            _dos->emplace_back(DedicatedObserver<_System, CHOSEN_STATE>());        
            initDOs<n-1, DOs>::init(_dos, _sys);
        }
    };
    template <typename DOs>
    struct initDOs<0, DOs>{
        static constexpr std::size_t CHOSEN_STATE{_System::n_outputs-1};
        static auto init(DOs* _dos, _System* _sys){
            _dos->emplace_back(DedicatedObserver<_System, CHOSEN_STATE>());
        }
    };
    template <int n, typename DOs>
    struct setSubjectDOs{
        static auto set(DOs* _dos, _System* _sys){
            //-- idk why with boost::variant the virtual method won't override
            ::jacl::system::detail::BaseSystemClient<typename _System::base_t>::setSubject(
                &boost::get<DedicatedObserver<_System, n>>((*_dos)[n]),
                ::jacl::system::detail::BaseSystemClient<typename _System::base_t>::ss(_sys));
            setSubjectDOs<n-1, DOs>::set(_dos, _sys);
        }
    };
    template <typename DOs>
    struct setSubjectDOs<0, DOs>{
        static auto set(DOs* _dos, _System* _sys){
            ::jacl::system::detail::BaseSystemClient<typename _System::base_t>::setSubject(
                &boost::get<DedicatedObserver<_System, 0>>((*_dos)[0]),
                ::jacl::system::detail::BaseSystemClient<typename _System::base_t>::ss(_sys));
        }
    };
    struct ifd_tag{};
    struct sifd_tag{};
    template <int n,typename DOs>
    struct setPoleDOs{
        static auto set(DOs* _dos,
                        const std::vector<typename arma::Col<
                            typename _System::scalar_t>::template fixed<_System::n_states-1>>& _pole, ifd_tag){
            boost::get<DedicatedObserver<_System, n>>( (*_dos)[n] ).setPole(_pole[n]);
            setPoleDOs<n-1,DOs>::set(_dos, _pole, ifd_tag());
        }
        static auto set(DOs* _dos,
                        const std::vector<typename arma::Col<
                            typename _System::scalar_t>::template fixed<_System::n_states-1>>& _pole, sifd_tag){
            boost::get<DedicatedObserver<_System, n>>( (*_dos)[n] ).setPole(_pole[n]);
        }
    };
    template <typename DOs>
    struct setPoleDOs<0, DOs>{
        static auto set(DOs* _dos,
                        const std::vector<typename arma::Col<
                            typename _System::scalar_t>::template fixed<_System::n_states-1>>& _pole, ifd_tag){
            boost::get<DedicatedObserver<_System, 0>>( (*_dos)[0] ).setPole(_pole[0]);
        }
        static auto set(DOs* _dos,
                        const std::vector<typename arma::Col<
                            typename _System::scalar_t>::template fixed<_System::n_states-1>>& _pole, sifd_tag){
            boost::get<DedicatedObserver<_System, 0>>( (*_dos)[0] ).setPole(_pole[0]);
        }
    };
    template <int n, typename DOs>
    struct getEst{
        static auto get(DOs* _dos,
                        const typename _System::input_t& _in,
                        const typename _System::output_t& _meas,
                        std::vector<typename _System::output_t>* _y_hat){
            static constexpr std::size_t IDX{_System::n_outputs-1 - n};
            _y_hat->emplace_back(boost::get<DedicatedObserver<_System, IDX>>( (*_dos)[IDX] ).convolve(_in, _meas));
            getEst<n-1, DOs>::get(_dos, _in, _meas, _y_hat);
        }
        static auto get(DOs* _dos,
                        const typename _System::input_t& _in,
                        const typename _System::output_t& _meas,
                        typename _System::output_t* _y_hat){
            *_y_hat = boost::get<DedicatedObserver<_System, n>>( (*_dos)[n] ).convolve(_in, _meas);
        }
    };
    template <typename DOs>
    struct getEst<0,DOs>{
        static auto get(DOs* _dos,
                        const typename _System::input_t& _in,
                        const typename _System::output_t& _meas,
                        std::vector<typename _System::output_t>* _y_hat){
            static constexpr std::size_t IDX{_System::n_outputs-1 - 0};
            _y_hat->emplace_back(boost::get<DedicatedObserver<_System, IDX>>( (*_dos)[IDX] ).convolve(_in, _meas));
        }
        static auto get(DOs* _dos,
                        const typename _System::input_t& _in,
                        const typename _System::output_t& _meas,
                        typename _System::output_t* _y_hat){
            *_y_hat = boost::get<DedicatedObserver<_System, 0>>( (*_dos)[0] ).convolve(_in, _meas);
        }
    };
    template <std::size_t num_sig, std::size_t ...sig_idx>
    struct plotSig{
        static auto idx()
            -> std::initializer_list<std::size_t>{            
            return plotSig<num_sig-1,
                        _System::n_states
                            + _System::n_inputs
                            + _System::n_outputs
                            - num_sig, sig_idx...>::idx();
        }
    };
    template <std::size_t ...sig_idx>
    struct plotSig<0, sig_idx...>{
        static auto idx()
            -> std::initializer_list<std::size_t>{            
            return {_System::n_states
                            + _System::n_inputs
                            + _System::n_outputs,
                            sig_idx...};
        }
    };
};

} } // namespace jacl::diagnosis