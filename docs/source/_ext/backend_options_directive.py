"""Sphinx directive rendering EDRIXS backend option registries."""

from docutils import nodes
from docutils.parsers.rst import Directive

from edrixs._backend_options import default_text
from edrixs.dense_backend.options import OPTIONS as DENSE_OPTIONS
from edrixs.fortran_backend.options import OPTIONS as FORTRAN_OPTIONS
from edrixs.petsc_backend.options import OPTIONS as PETSC_OPTIONS
from edrixs.scipy_backend.options import OPTIONS as SCIPY_OPTIONS


REGISTRIES = {
    'dense': DENSE_OPTIONS,
    'fortran': FORTRAN_OPTIONS,
    'petsc': PETSC_OPTIONS,
    'scipy': SCIPY_OPTIONS,
}


def _entry(text, *, literal=False):
    entry = nodes.entry()
    paragraph = nodes.paragraph()
    paragraph += nodes.literal(text=text) if literal else nodes.Text(text)
    entry += paragraph
    return entry


def _row(values, *, header=False):
    row = nodes.row()
    for index, value in enumerate(values):
        row += _entry(value, literal=not header and index < 2)
    return row


class BackendOptionsDirective(Directive):
    """Render the option table for ``BACKEND OPERATION``."""

    required_arguments = 2

    def run(self):
        backend, operation = self.arguments
        try:
            specifications = REGISTRIES[backend][operation]
        except KeyError as exc:
            raise self.error(
                "unknown backend option table: {} {}".format(
                    backend, operation
                )
            ) from exc

        if not specifications:
            return [nodes.paragraph(text='No backend-specific options.')]

        table = nodes.table()
        tgroup = nodes.tgroup(cols=3)
        table += tgroup
        for width in (20, 24, 56):
            tgroup += nodes.colspec(colwidth=width)

        head = nodes.thead()
        head += _row(('Option', 'Default', 'Description'), header=True)
        tgroup += head

        body = nodes.tbody()
        for name, specification in specifications.items():
            body += _row((
                name,
                default_text(specification),
                specification.description,
            ))
        tgroup += body
        return [table]


def setup(app):
    app.add_directive('backend-options', BackendOptionsDirective)
    return {
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }
