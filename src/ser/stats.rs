use crate::Result;
use csv::WriterBuilder;

#[derive(serde::Serialize)]
#[serde(rename_all = "snake_case")]
pub enum StatsStage {
    Resample,
    ResampleRegion,
    Seq,
    Climb,
}

#[derive(serde::Serialize)]
pub struct StatsRecord {
    pub stage: StatsStage,
    pub iter: usize,
    pub prior: f64,
}

pub struct StatsWriter {
    _writer: csv::Writer<std::fs::File>,
}

impl StatsWriter {
    pub fn from_path<P: AsRef<std::path::Path>>(path: P) -> Result<Self> {
        let file = std::fs::File::create(path)?;
        let writer = WriterBuilder::new().delimiter(b'\t').from_writer(file);
        Ok(Self { _writer: writer })
    }
    pub fn serialize(&mut self, stats: &StatsRecord) -> Result<()> {
        Ok(self._writer.serialize(stats)?)
    }
}
