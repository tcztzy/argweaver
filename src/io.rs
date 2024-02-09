#[cfg(feature = "extension-module")]
use pyo3::{prelude::*, types::PyDict};
pub struct Stdout;

impl std::io::Write for Stdout {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        let s = std::str::from_utf8(buf).unwrap();
        #[cfg(feature = "extension-module")]
        Python::with_gil(|py| {
            let locals = PyDict::new(py);
            locals.set_item("s", s).unwrap();
            py.run("import sys;sys.stdout.write(s)", None, Some(&locals))
                .unwrap();
        });
        #[cfg(not(feature = "extension-module"))]
        write!(std::io::stdout(), "{}", s).unwrap();
        Ok(buf.len())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        #[cfg(feature = "extension-module")]
        {
            Python::with_gil(|py| {
                py.run("import sys;sys.stdout.flush()", None, None).unwrap();
            });
            Ok(())
        }
        #[cfg(not(feature = "extension-module"))]
        std::io::stdout().flush()
    }
}
