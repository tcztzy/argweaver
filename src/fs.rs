use flate2::read::GzDecoder;
use std::fs::File;
use std::io::Read;
use std::path::Path;

use crate::Result;

pub fn read_to_string(path: &Path) -> Result<String> {
    let mut file = File::open(path)?;
    let mut contents = String::new();
    if path.extension().unwrap() == "gz" {
        let mut gz = GzDecoder::new(file);
        gz.read_to_string(&mut contents)?;
    } else {
        file.read_to_string(&mut contents)?;
    }
    Ok(contents)
}
