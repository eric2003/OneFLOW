# time_integration/factory.py

from factories.base_factory import BaseFactory

class TimeIntegratorFactory:
    @staticmethod
    def create(cfd) -> 'TimeIntegrator':
        rk_order = cfd.config.rk_order
        integrator_name = f'rk{rk_order}'
        return BaseFactory.create_component('integrator', integrator_name, cfd)