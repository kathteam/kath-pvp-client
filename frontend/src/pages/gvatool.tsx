import { JSX } from "react"
import { useNavigate } from "react-router-dom"

export default function GVATool(): JSX.Element {
  const navigate = useNavigate()
  
  return (
    <div>
      <h1>GVATool Page</h1>
      <p>This is the gene variation analysis tool page of our application.</p>
      <button onClick={() => navigate('/index.html')}>
        Back to Dashboard
      </button>
    </div>
  )
}
