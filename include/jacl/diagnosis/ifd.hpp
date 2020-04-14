#pragma once

#include <boost/variant.hpp>

#include <jacl/traits.hpp>
#include <jacl/system/continuous_system.hpp>
#include <jacl/system/discrete_system.hpp>
#include <jacl/pole_placement.hpp>
#include <jacl/lti_common.hpp>

namespace jacl{ namespace diagnosis{

template <typename _System>
class IFD;

namespace detail{
    template <typename _System>
    class BaseSystemClient{
    protected:
        static auto ss(_System* _sys) -> typename _System::StateSpace*{
            if(_sys)
                return _sys->ss_;
            return nullptr;
        }
        template <typename __System, typename _StateSpace>
        static auto setSubject(__System* _sys, _StateSpace* _ss) -> void{
            _sys->setSubject(_ss);
        }
    };
} // namespace detail

template <typename _System>
class IFD:public detail::BaseSystemClient<typename _System::base>{
public:    
    IFD(_System* _sys, std::initializer_list<double> _threshold)
        : sys_(_sys){
        // , dos_(_System::n_outputs, DedicatedObserver<_System, 0>(sys_)){     
        
        auto _th_it = _threshold.begin();
        for(int i(0); i < threshold_.size() && _th_it != _threshold.end(); i++){
            threshold_[i] = *_th_it++;
        }
     
        initDOs<_System::n_outputs-1, decltype(dos_)>::init(&dos_, sys_);
        std::cout << "DOS Size : " << dos_.size() << std::endl;
        setSubjectDOs<_System::n_outputs-1, decltype(dos_)>::set(&dos_, sys_);
        printTest<_System::n_outputs-1, decltype(dos_)>::pr(&dos_);
    }
    ~IFD(){}
    virtual auto init(std::initializer_list<arma::vec > _poles) -> void{
        std::vector<arma::vec> poles(_poles);
        setPoleDOs<_System::n_outputs-1, decltype(dos_)>::set(&dos_, poles);
        // int i(0);
        // for(auto o:dos_){
        //     boost::apply_visitor([poles,i,o](decltype(o) a){a.setPole(poles[i]);},o);
        //     ++i;
        // }
        // for(int i(0); i < _System::n_states; i++){            
        //     dos_[i]->setPole(poles[i]);            
        // }
    }
    virtual auto detect(const arma::vec& _in, const arma::vec& _out)
        -> std::array<std::pair<arma::vec, bool>, _System::n_outputs>{
        std::array<std::pair<arma::vec, bool>, _System::n_outputs> res;
        std::vector<arma::vec> y_hat;
        arma::vec delta_y;
        getEst<_System::n_outputs-1, decltype(dos_)>::get(&dos_, _in, _out, &y_hat);
        for(int i(0); i < _System::n_outputs; i++){
            delta_y = _out - y_hat[i];
            auto err(1.);
            for(int j(0); j < _System::n_outputs; j++){
                if(j != i)
                    err *= delta_y(j);
            }
            std::get<0>(res[i]) = y_hat[i];
            std::get<1>(res[i]) = err > threshold_[i];
        }
        return res;
    }

protected:
    template <typename __System = _System, std::size_t chosen_state=0>
    class DedicatedObserver:public ::jacl::system::BaseSystem<typename __System::StateSpace>{
    public:
        static constexpr double CHOSEN_STATE{chosen_state};
        DedicatedObserver(__System* _sys=nullptr)
            : ::jacl::system::BaseSystem<typename __System::StateSpace>(
                    detail::BaseSystemClient<typename __System::base>::ss(_sys))
            , z_( arma::zeros<arma::vec>(__System::n_outputs-1, 1) )
            , prev_z_(z_){
            std::cout << "DOS-" << CHOSEN_STATE << " created !!" << std::endl;
        }
        // DedicatedObserer()
        //     : z_( arma::zeros<arma::vec>(__System::n_outputs-1, 1) )
        //     , prev_z_(z_){
        // }
        ~DedicatedObserver(){}
         //-- overload
        auto convolve(const arma::vec& _in, const arma::vec& _meas) -> arma::vec{
            meas_ = _meas;
            arma::vec est( convolve(_in) );
            prev_meas_ = meas_;
            return est;
        }    
        void setPole(const arma::vec& _poles){
            poles_ = _poles;
            updateVar();
        }
        void printTest(){
            std::cout << "DOS-" << CHOSEN_STATE << " masuk !!" << std::endl;
            std::cout << this->ss_ << std::endl;
            // std::cout << &(*this->ss_) << std::endl;
        }
    protected:
        auto convolve(const arma::vec& _in) -> arma::vec override{
            setIn(_in);
            this->state_ = dstate();        
            this->out_ = output();
            this->prev_state_ = this->state_;
            return this->out_;
        }
        auto setIn(const arma::vec& _in) -> void{
            this->in_ = _in;
        }
        template <typename T = __System>
        auto dstate(const arma::mat& term1,
                    const arma::mat& term2,
                    const arma::mat& term3,
                    const arma::vec& _yI)
            -> typename std::enable_if<
                ::jacl::traits::is_continuous_system<T>::value, arma::vec>::type{
            //-- not implemented yet
            return arma::vec(T::n_outputs, 1, arma::fill::zeros);
        }
        template <typename T = __System>
        auto dstate(const arma::mat& _term1,
                    const arma::mat& _term2,
                    const arma::mat& _term3,
                    const arma::vec& _yI)
            -> typename std::enable_if<
                ::jacl::traits::is_discrete_system<T>::value, arma::vec>::type{
            return _term1*prev_z_ + _term1*L_*_yI
                    + _term2*_yI + _term3*this->in_; 
        }
        auto dstate() -> arma::vec override{
            arma::mat term1( A22_ - L_*A12_ );
            arma::mat term2( A21_ - L_*A11_ );
            arma::mat term3( B2_ - L_*B1_ );
            arma::vec yI( meas_(CHOSEN_STATE) - this->ss_->D().row(CHOSEN_STATE) * this->in_ );
            z_ = dstate(term1, term2, term3, yI);
            q_hat_ = arma::join_vert(yI, prev_z_ + L_*yI);
            prev_z_ = z_;
            return  arma::inv(M_) * q_hat_;
        }
        auto output() -> arma::vec override{
            return this->ss_->C() * this->prev_state_ + this->ss_->D() * this->in_;
        }
        auto updateVar() -> void override{
            if(!An_.is_empty() && !Bn_.is_empty()){
                transform();
                An_.print("An : ");
                Bn_.print("Bn : ");
                A11_ = A11(An_);
                A12_ = A12(An_);
                A21_ = A21(An_);
                A22_ = A22(An_);
                B1_ = B1(Bn_);
                B2_ = B2(Bn_);
                jacl::LinearStateSpace<_System::n_states-1,1,1> ss(A22_, arma::zeros(_System::n_states-1,1),
                                                                    A12_, arma::zeros(1,1));
                jacl::pole_placement::KautskyNichols(&ss, poles_, &L_, jacl::pole_placement::PolePlacementType::Observer);

                L_.print("Gain : ");
                // if(jacl::common::observable(ss.A(), ss.C()))
                //     std::cout << "SS Observable !" << std::endl;
                // arma::cx_mat eigvec;
                // arma::cx_vec eigval;
                // arma::eig_gen(eigval, eigvec, ss.A() - gain*ss.C());
                // eigval.print("EigVal : ");
            }
        }
        auto transform() -> void{
            std::cout << "ASD1" << std::endl;
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
            std::cout << "ASD2" << std::endl;     
        }
        inline auto A11(const arma::mat& _A) -> arma::mat{
            return _A.submat(0,0,0,0);
        }
        inline auto A12(const arma::mat& _A) -> arma::mat{
            return _A.submat(0,1,0,__System::n_states-1);
        }
        inline auto A21(const arma::mat& _A) -> arma::mat{
            return _A.submat(1,0,__System::n_states-1,0);
        }
        inline auto A22(const arma::mat& _A) -> arma::mat{
            return _A.submat(1,1, __System::n_states-1 ,__System::n_states-1);
        }
        inline auto B1(const arma::mat& _B) -> arma::mat{
            return _B.head_rows(1);
        }
        inline auto B2(const arma::mat& _B) -> arma::mat{
            return _B.tail_rows(__System::n_states-1);
        }
    protected:
        arma::mat An_;
        arma::mat Bn_;
        arma::mat A11_;
        arma::mat A12_;
        arma::mat A21_;
        arma::mat A22_;
        arma::mat B1_;
        arma::mat B2_;
        arma::mat M_;        
        arma::mat L_;
        arma::vec z_;
        arma::vec prev_z_;
        arma::vec q_hat_;
        arma::vec meas_;
        arma::vec prev_meas_;
        arma::vec poles_;
    };

protected:
    std::array<double, _System::n_outputs> threshold_;
    _System* sys_;

    template <int n, typename T>
    struct DOVariant;
    template <int n, typename ...Args>
    struct DOVariant<n, boost::variant<Args...> >{
        using type = typename DOVariant<n-1, boost::variant<Args..., DedicatedObserver<_System,n-1> > >::type;
    };
    template <typename ...Args>
    struct DOVariant<1, boost::variant<Args...> >{
        using type = typename boost::variant<Args..., DedicatedObserver<_System,0> >;
    };
    // using DOS = std::vector<DOVariant<_System::n_outputs-1, boost::variant<DedicatedObserver<_System, _System::n_outputs-1> > >::t >;
    using DOS = std::vector<boost::variant<DedicatedObserver<_System,0>, DedicatedObserver<_System,1>, DedicatedObserver<_System,2> > >;
    // DOS dos_{DedicatedObserver<_System,0>(sys_),DedicatedObserver<_System,1>(sys_),DedicatedObserver<_System,2>(sys_)};
    DOS dos_;
    
    template <int n, typename DOs>
    struct initDOs{
        static constexpr std::size_t CHOSEN_STATE{(_System::n_outputs-1) - n};
        static auto init(DOs* _dos, _System* _sys) -> void{
            _dos->emplace_back(DedicatedObserver<_System, CHOSEN_STATE>());
            //-- idk why with boost::variant the virtual method won't override
            // detail::BaseSystemClient<typename _System::base>::setSubject(
            //     &boost::get<DedicatedObserver<_System, CHOSEN_STATE>>(_dos->back()), detail::BaseSystemClient<typename _System::base>::ss(_sys));
            // boost::get<DedicatedObserver<_System, n>>(_dos->back()).setSubject(detail::BaseSystemClient<typename _System::base>::ss(_sys));
            initDOs<n-1, DOs>::init(_dos, _sys);
        }
    };
    template <typename DOs>
    struct initDOs<0, DOs>{
        static constexpr std::size_t CHOSEN_STATE{_System::n_outputs-1};
        static auto init(DOs* _dos, _System* _sys) -> void{
            _dos->emplace_back(DedicatedObserver<_System, CHOSEN_STATE>());
            // detail::BaseSystemClient<typename _System::base>::setSubject(
            //     &boost::get<DedicatedObserver<_System, CHOSEN_STATE>>(_dos->back()), detail::BaseSystemClient<typename _System::base>::ss(_sys));
            // boost::get<DedicatedObserver<_System, 0>>(_dos->back()).setSubject(detail::BaseSystemClient<typename _System::base>::ss(_sys));
        }
    };
    template <int n, typename DOs>
    struct setSubjectDOs{
        static auto set(DOs* _dos, _System* _sys) -> void{
            detail::BaseSystemClient<typename _System::base>::setSubject(
                &boost::get<DedicatedObserver<_System, n>>((*_dos)[n]), detail::BaseSystemClient<typename _System::base>::ss(_sys));
            setSubjectDOs<n-1, DOs>::set(_dos, _sys);
        }
    };
    template <typename DOs>
    struct setSubjectDOs<0, DOs>{
        static auto set(DOs* _dos, _System* _sys) -> void{
            detail::BaseSystemClient<typename _System::base>::setSubject(
                &boost::get<DedicatedObserver<_System, 0>>((*_dos)[0]), detail::BaseSystemClient<typename _System::base>::ss(_sys));
        }
    };
    template <int n,typename DOs>
    struct setPoleDOs{
        static auto set(DOs* _dos, const std::vector<arma::vec>& _pole) -> void{
            boost::get<DedicatedObserver<_System, n>>( (*_dos)[n] ).setPole(_pole[n]);
            setPoleDOs<n-1,DOs>::set(_dos, _pole);
        }
    };
    template <typename DOs>
    struct setPoleDOs<0, DOs>{
        static auto set(DOs* _dos, const std::vector<arma::vec>& _pole) -> void{
            boost::get<DedicatedObserver<_System, 0>>( (*_dos)[0] ).setPole(_pole[0]);
        }
    };
    template <int n, typename DOs>
    struct getEst{
        static auto get(DOs* _dos, const arma::vec& _in, const arma::vec& _meas, std::vector<arma::vec>* _y_hat) -> void{
            _y_hat->emplace_back(boost::get<DedicatedObserver<_System, n>>( (*_dos)[n] ).convolve(_in, _meas));
            getEst<n-1, DOs>::get(_dos, _in, _meas, _y_hat);
        }
    };
    template <typename DOs>
    struct getEst<0,DOs>{
        static auto get(DOs* _dos, const arma::vec& _in, const arma::vec& _meas, std::vector<arma::vec>* _y_hat) -> void{
            _y_hat->emplace_back(boost::get<DedicatedObserver<_System, 0>>( (*_dos)[0] ).convolve(_in, _meas));
        }
    };
    template <int n,typename DOs>
    struct printTest{
        static auto pr(DOs* _dos) -> void{
            boost::get<DedicatedObserver<_System, n>>( (*_dos)[n] ).printTest();
            printTest<n-1,DOs>::pr(_dos);
        }
    };
    template <typename DOs>
    struct printTest<0, DOs>{
        static auto pr(DOs* _dos) -> void{
            boost::get<DedicatedObserver<_System, 0>>( (*_dos)[0] ).printTest();
        }
    };       
};

} } // namespace jacl::diagnosis