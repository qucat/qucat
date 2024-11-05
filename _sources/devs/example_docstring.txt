def zpf(self, mode, quantity, **kwargs):
    r'''Returns contribution of a mode to the zero-point fluctuations of a quantity for this component.

    The quantity can be current (in units of Ampere), 
    voltage (in Volts), 
    charge (in electron charge), 
    or flux (in units of the reduced flux quantum, :math:`\hbar/2e`).

    Parameters
    ----------
    mode:           integer
                    Determine what mode to consider, where 0 designates
                    the lowest frequency mode, and the others
                    are arranged in order of increasing frequency
    quantity:       string
                    One of 'current', 'flux', 'charge', 'voltage'
    kwargs:     
                Values for un-specified circuit components, 
                ex: ``L=1e-9``.

    Returns
    -------
    float
        contribution of the ``mode`` to the zero-point fluctuations of the ``quantity``

    Notes
    -----
    This quantity is calculated by multiplying the
    voltage transfer function :math:`T_{rc}` (between a reference component :math:`r`
    and the annotated component  :math:`c` ), with
    :math:`X_{zpf,m,r}`, the zero-point fluctuations of :math:`\hat{X}` at the reference component.
    
    Note that resistors make the transfer function :math:`T_{rc}`, and hence this quantity, complex.

    For more detail on the underlying theory, see https://arxiv.org/pdf/1908.10342.pdf.
    ''''''