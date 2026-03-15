#ifndef LIMITERCONTROLLER_H
#define LIMITERCONTROLLER_H
namespace limiter
{
    void limiter_1InnerStep();

    void limiter_1OutterStep();

    void limitRho_PositivityPreserving();

    namespace IOLimiter {
        void readSelectedLimiters();
    }
}
#endif // LIMITERCONTROLLER_H
