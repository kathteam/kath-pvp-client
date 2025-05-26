import { Box, Button, Modal, Typography } from "@mui/material";
import { JSX, useState } from "react";

// Import all GIF files from the Media directory
import extractButtonGif from "../../Media/ExtractButton.gif";
import extractStringGif from "../../Media/ExtractString.gif";
import geneDownloadGif from "../../Media/GeneDownload.gif";
import downloadFastaGif from "../../Media/DownloadFasta.gif";

// Map of media names to their imported paths
const mediaMap: { [key: string]: string } = {
	"ExtractButton.gif": extractButtonGif,
	"ExtractString.gif": extractStringGif,
	"GeneDownloading.gif": geneDownloadGif,
	"DownloadFasta.gif": downloadFastaGif,
};

interface HowToProps {
	media?: string;
	title?: string;
	description?: string;
}

export default function HowTo({
	media = "ExtractButton.gif",
	title = "How to use",
	description = "Click the button to see the how to use",
}: HowToProps): JSX.Element {
	const [open, setOpen] = useState(false);

	const handleOpen = () => setOpen(true);
	const handleClose = () => setOpen(false);

	const isGif = media.toLowerCase().endsWith(".gif");
	const mediaUrl = mediaMap[media] || "";

	return (
		<>
			<Button variant="outlined" color="secondary" onClick={handleOpen}>
				How to use
			</Button>
			<Modal open={open} onClose={handleClose}>
				<Box
					sx={{
						position: "absolute",
						top: "50%",
						left: "50%",
						transform: "translate(-50%, -50%)",
						width: "40vw",
						bgcolor: "background.paper",
						boxShadow: 24,
						p: 4,
						display: "flex",
						flexDirection: "column",
						alignItems: "center",
						borderRadius: 2,
					}}
				>
					{isGif ? (
						<img
							src={mediaUrl}
							alt={title}
							style={{ width: "100%", height: "auto" }}
						/>
					) : (
						<video
							src={mediaUrl}
							autoPlay
							loop
							muted
							style={{ width: "100%", height: "auto" }}
						/>
					)}
					<Typography variant="h5" component="h2" sx={{ mb: 3, mt: 2 }}>
						{title}
					</Typography>
					<Typography variant="body1" component="p" sx={{ mb: 3 }}>
						{description}
					</Typography>
					<Button
						variant="contained"
						color="primary"
						onClick={handleClose}
						sx={{ mt: 2 }}
					>
						Close
					</Button>
				</Box>
			</Modal>
		</>
	);
}
