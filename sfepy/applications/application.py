from sfepy.base.base import Struct, output, insert_as_static_method

class Application(Struct):
    """
    Base class for applications.

    Subclasses should implement: __init__(), call().

    Automates parametric studies, see parametrize().
    """
    def __init__(self, conf, options, output_prefix, **kwargs):
        Struct.__init__(self,
                        conf=conf,
                        options=options,
                        output_prefix=output_prefix)
        output.prefix = self.output_prefix
        self.restore()

    def setup_options(self):
        pass

    def __call__(self, **kwargs):
        """
        This is either call_basic() or call_parametrized().
        """
        pass

    def call_basic(self, **kwargs):
        return self.call(**kwargs)

    def call_parametrized(self, **kwargs):
        generator = self.parametric_hook(self.problem)
        for aux in generator:
            if isinstance(aux, tuple) and (len(aux) == 2):
                problem, container = aux
                mode = 'coroutine'
            else:
                problem = aux
                mode = 'simple'
            self.problem = problem

            generator_prefix = output.prefix
            output.prefix = self.output_prefix # Restore default.

            """Application options have to be re-processed here as they can
            change in the parametric hook."""
            self.setup_options()
            out = self.call(**kwargs)

            output.prefix = generator_prefix

            if mode == 'coroutine':
                # Pass application output to the generator.
                container.append(out)
                next(generator)

    def restore(self):
        """
        Remove parametric_hook, restore __call__() to call_basic().
        """
        self.parametric_hook = None
        insert_as_static_method(self.__class__, '__call__',
                                self.call_basic)

    def parametrize(self, parametric_hook):
        """
        Add parametric_hook, set __call__() to call_parametrized().
        """
        if parametric_hook is None: return

        self.parametric_hook = parametric_hook
        insert_as_static_method(self.__class__, '__call__',
                                self.call_parametrized)
